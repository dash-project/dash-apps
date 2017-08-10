#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cstdio>

#include <libdash.h>
#include <dash/experimental/HaloMatrixWrapper.h>

#include "minimonitoring.h"

//#define WITHPNGOUTPUT 1
#ifdef WITHPNGOUTPUT
#include <png++/png.hpp>

uint32_t filenumber= 0;
const std::string filename("image");
const std::string suffix("image");

#endif /* WITHPNGOUTPUT */


/* for MiniMonT */
std::vector<MiniMonT>* MiniMonT::tape;
std::chrono::time_point<std::chrono::high_resolution_clock> MiniMonT::inittime;


using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::vector;

using TeamSpecT = dash::TeamSpec<2>;
using MatrixT = dash::NArray<double,2>;
using StencilT = dash::experimental::Stencil<2>;
using StencilSpecT = dash::experimental::StencilSpec<2,4>;
using CycleSpecT = dash::experimental::CycleSpec<2>;
using HaloMatrixWrapperT = dash::experimental::HaloMatrixWrapper<MatrixT,StencilSpecT>;


const StencilSpecT stencil_spec({
    StencilT(-1, 0), StencilT( 1, 0),
    StencilT( 0,-1), StencilT( 0, 1)});

const CycleSpecT cycle_spec( dash::experimental::Cycle::FIXED, dash::experimental::Cycle::FIXED );

struct Level {

    MatrixT grid;
    HaloMatrixWrapperT halo;

    Level( size_t h, size_t w, dash::Team& team, TeamSpecT teamspec ) :
        grid( dash::SizeSpec<2>( h, w ),
            dash::DistributionSpec<2>( dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
        halo( grid, stencil_spec, cycle_spec ) {}
    Level() = delete;

    ~Level() {}
};

#ifdef WITHPNGOUTPUT
void topng( const MatrixT& grid, double minvalue= 0.0, double maxvalue= 100.0 ) {

    if ( 0 != dash::myid() ) return;

    uint32_t minsize= 512;
    uint32_t f= 1;
    size_t w= grid.extent(1);
    size_t h= grid.extent(0);

    if ( w < minsize || h < minsize ) {

        f= ( w <= h ) ? w : h;
        f= 1+minsize/f;
    }

    png::image< png::rgb_pixel > image( w*f, h*f );
    for ( png::uint_32 y = 0; y < h; y++ ) {
        for ( png::uint_32 x = 0; x < w; x++ ) {

            double temp= grid[y][x];
            double r= (temp-minvalue)*255.0/(maxvalue-minvalue);
            double g= 0;
            double b= 255.0 - r;

            r= ( r > 255.0 ) ? 255.0 : ( r < 0.0 ) ? 0.0 : r ;
            g= ( g > 255.0 ) ? 255.0 : ( g < 0.0 ) ? 0.0 : g ;
            b= ( b > 255.0 ) ? 255.0 : ( b < 0.0 ) ? 0.0 : b ;

            for ( png::uint_32 yy = 0; yy < f; yy++ ) {
                for ( png::uint_32 xx = 0; xx < f; xx++ ) {
                    image[y*f+yy][x*f+xx] = png::rgb_pixel( (uint8_t) r, (uint8_t) g, (uint8_t) b );
                }
            }
        }
    }

    char buf[10];
    std::sprintf( buf, "%03u", filenumber++ );

    image.write( std::string( filename + std::string(buf) + std::string(".png") ).c_str() );
}
#endif /* WITHPNGOUTPUT */

void writeToPNG(const MatrixT& grid) {

#ifdef WITHPNGOUTPUT
    grid.barrier();

    if( 0 == dash::myid() )
        topng( grid, 0.0, 10.0 );

    grid.barrier();
#endif /* WITHPNGOUTPUT */
}

void sanitycheck( const MatrixT& grid  ) {

    /* check if the sum of the local extents of the matrix blocks sum up to
    the global extents, abort otherwise */
    dash::Array<size_t> sums( dash::Team::All().size(), dash::BLOCKED );
    sums.local[0]= grid.local.extent(0) * grid.local.extent(1);
    sums.barrier();

    if ( 0 != dash::myid() ) {

        dash::transform<size_t>(
            sums.lbegin(), sums.lend(), // first source
            sums.begin(), // second source
            sums.begin(), // destination
            dash::plus<size_t>() );
    }
    sums.barrier();

    if ( ( (size_t) sums[0] ) != grid.extent(0) * grid.extent(1) ) {

        if ( 0 == dash::myid() ) {
            cout << "ERROR: size mismatch: global size is " <<
                grid.extent(0) * grid.extent(1) << " == " << grid.extent(0) << "x" << grid.extent(1) <<
                " but sum of local sizes is " << (size_t) sums[0] << std::flush << endl;
        }

        dash::finalize();
        exit(0);
    }
}


void initgrid( MatrixT& grid ) {

    /* init to 0 first! */
    dash::fill( grid.begin(), grid.end(), 0.0 );

    if ( 0 == dash::myid() ) {

        /* do initialization only from unit 0 with global accesses.
        This is not the fastest way and a parallel routine would not be too
        complicated but it should demonstrate the possibility to trade
        convenience against performance. */

        size_t w= grid.extent(1);
        size_t h= grid.extent(0);

        for( size_t i = 0; i < h; ++i ) {
            for( size_t k = 0; k < w; ++k ) {
                grid[i][k] = 5.0;
            }
        }

        for( size_t i = h/5; i < h*3/4; ++i ) {
            for( size_t k = w*1/5; k < w*1/3; ++k ) {
                grid[i][k] = 3.0;
            }
        }

    }

    grid.barrier();
}


void initboundary( Level& level ) {

    double gw= level.grid.extent(1);
    double gh= level.grid.extent(0);

    /* with this way of setting the boundary conditions on every level separatel (in
    contrast to scaling it down step by step) one needs to make sure that the boundary
    values are continuous (kind of). Otherwise the coarsened version of the boundary
    values will have jumps in different places (instead of jumps being smoothed out
    by coarsening the boundary).
    Well, we should really think about another way to init the boundary values! */

    level.halo.setFixedHalos( [gh,gw]( const std::array<dash::default_index_t,2>& coords ) {

        dash::default_index_t y= coords[0];
        dash::default_index_t x= coords[1];
        double ret= 0;
        if ( -1 == y ) {

            /* upper border */
            if ( x / gw <= 0.5 ) ret= 10.0 - 10.0*(x+0.5)/gw/0.5;
            /* the +0.5 is important to make the next level's boundary the same
            as the previous one scaled down. Only then does the residual behave
            like expected! Same for the other 3 cases following. */

        } else if ( gh == y ) {

            /* lower border */
            if ( x / gw > 0.4 ) ret= 10.0 - 10.0*(gw-(x+0.5))/gw/0.6;

        } else if ( -1 == x ) {

            /* left border */
            if ( y / gh <=0.5 ) ret= 10.0 - 10.0*(y+0.5)/gh/0.5;

        } else if ( gw == x ) {

            /* right border */
            if ( y / gh > 0.8 ) ret= 10.0 - 10.0*(gh-(y+0.5))/gh/0.2;
        }

        /*
        cerr << "set " << coords[0] << "/" << gh << "," << coords[1] << "/" << gw << "== " <<
            coords[0]/gh << "," << coords[1]/gw << ":" << ret << endl;
        */

        return ret;
    });

}


void markunits( MatrixT& grid ) {

    /* Mark unit bordery by setting the first local rows and columns */

    size_t w= grid.local.extent(1);
    size_t h= grid.local.extent(0);

    for( size_t i = 0; i < h; ++i ) {
        grid.local[i][0] = 3.0;
    }

    for( size_t k = 0; k < w; ++k ) {
        grid.local[0][k] = 3.0;
    }

}


void markboundary( Level& level ) {

    auto local= level.grid.lbegin();

    auto bend = level.halo.bend();
    for( auto it = level.halo.bbegin(); it != bend; ++it ) {

        double dtheta = ( it.valueAt(2) + it.valueAt(3) +
            it.valueAt(0) + it.valueAt(1) ) / 4.0 - *it;
        local[ it.lpos() ] += dtheta;
    }
}


void scaledownboundary( const MatrixT& fine, MatrixT& coarse ) {

    assert( (coarse.extent(1)-2)*2 +2 == fine.extent(1) );
    assert( (coarse.extent(0)-2)*2 +2 == fine.extent(0) );

    size_t wc= coarse.local.extent(1);
    size_t hc= coarse.local.extent(0);
    size_t wf= fine.local.extent(1);
    size_t hf= fine.local.extent(0);
    std::array< long int, 2 > cornerc= coarse.pattern().global( {0,0} );

    size_t startx= ( 0 == cornerc[1] ) ? 1 : 0;
    if ( 0 == cornerc[0] ) {

        /* if boundary in front */
        for ( size_t x= 0; x < wc-1; x++ ) {
            coarse.local[0][startx+x]= 0.5 * (double) ( fine.local[0][startx+2*x+0] + fine.local[0][startx+2*x+1] );
        }

    } else {

        /* else boundary in the end -- this is only true for the restructed case with 4 units */
        for ( size_t x= 0; x < wc-1; x++ ) {
            coarse.local[hc-1][startx+x]= 0.5 * ( fine.local[hf-1][startx+2*x+0] + fine.local[hf-1][startx+2*x+1] );
        }
    }

    size_t starty= ( 0 == cornerc[0] ) ? 1 : 0;
    if ( 0 == cornerc[1] ) {

        /* if boundary in front */
        for ( size_t y= 0; y < hc-1; y++ ) {
            coarse.local[starty+y][0]= 0.5 * ( fine.local[starty+2*y+0][0] + fine.local[starty+2*y+1][0] );
        }

    } else {

        /* else boundary in the end -- this is only true for the restructed case with 4 units */
        for ( size_t y= 0; y < hc-1; y++ ) {
            coarse.local[starty+y][wc-1]= 0.5 * ( fine.local[starty+2*y+0][wf-1] + fine.local[starty+2*y+1][wf-1] );
        }
    }
}


void scaledown( Level& fine, Level& coarse ) {

    uint32_t param= fine.grid.extent(0);
    MiniMonT::MiniMonRecord( 0, "scaledown", param );


    assert( coarse.grid.extent(1) * 2 == fine.grid.extent(1) );
    assert( coarse.grid.extent(0) * 2 == fine.grid.extent(0) );

    const std::array< long unsigned int, 2 > extentc= coarse.grid.pattern().local_extents();
    const std::array< long signed int, 2 >   cornerc= coarse.grid.pattern().global( {0,0} );
    const std::array< long unsigned int, 2 > extentf= fine.grid.pattern().local_extents();
    const std::array< long signed int, 2 >   cornerf= fine.grid.pattern().global( {0,0} );

    assert( cornerc[0] * 2 == cornerf[0] );
    assert( cornerc[1] * 2 == cornerf[1] );

    assert( 0 == cornerc[0] %2 );
    assert( 0 == cornerc[1] %2 );

    assert( extentc[0] * 2 == extentf[0] );
    assert( extentc[1] * 2 == extentf[1] );

    for ( size_t y= 0; y < extentc[0]; y++ ) {
        for ( size_t x= 0; x < extentc[1]; x++ ) {

            coarse.grid.local[y][x]= 0.25 * (
                fine.grid.local[2*y  ][2*x  ] +
                fine.grid.local[2*y  ][2*x+1] +
                fine.grid.local[2*y+1][2*x  ] +
                fine.grid.local[2*y+1][2*x+1] );
        }
    }

    MiniMonT::MiniMonRecord( 1, "scaledown", param );
}


void scaleup( Level& coarse, Level& fine ) {

    uint32_t param= coarse.grid.extent(0);
    MiniMonT::MiniMonRecord( 0, "scaleup", param );

    assert( coarse.grid.extent(1) * 2 == fine.grid.extent(1) );
    assert( coarse.grid.extent(0) * 2 == fine.grid.extent(0) );

    const std::array< long unsigned int, 2 > extentc= coarse.grid.pattern().local_extents();
    const std::array< long signed int, 2 >   cornerc= coarse.grid.pattern().global( {0,0} );
    const std::array< long unsigned int, 2 > extentf= fine.grid.pattern().local_extents();
    const std::array< long signed int, 2 >   cornerf= fine.grid.pattern().global( {0,0} );

    assert( cornerc[0] * 2 == cornerf[0] );
    assert( cornerc[1] * 2 == cornerf[1] );

    assert( 0 == cornerc[0] %2 );
    assert( 0 == cornerc[1] %2 );

    assert( extentc[0] * 2 == extentf[0] );
    assert( extentc[1] * 2 == extentf[1] );

    for ( size_t y= 0; y < extentc[0] ; y++ ) {
        for ( size_t x= 0; x < extentc[1]; x++ ) {

            double t= coarse.grid.local[y][x];
            fine.grid.local[2*y  ][2*x  ]= t;
            fine.grid.local[2*y  ][2*x+1]= t;
            fine.grid.local[2*y+1][2*x  ]= t;
            fine.grid.local[2*y+1][2*x+1]= t;
        }
    }

    MiniMonT::MiniMonRecord( 1, "scaleup", param );
}


double smoothen( Level& level ) {

    uint32_t param= level.grid.extent(0);
    MiniMonT::MiniMonRecord( 0, "smoothen", param );

    double res= 0.0;

    /// relaxation coeff.
    const double c= 1.0;

    /// begin pointer of local block, needed because halo border iterator is read-only
    auto gridlocalbegin= level.grid.lbegin();

    // async halo update
    level.halo.updateHalosAsync();

    // update inner

    MiniMonT::MiniMonRecord( 0, "smooth_inner", param );

    size_t lw= level.grid.local.extent(1);
    size_t lh= level.grid.local.extent(0);

    /* the start value for both, the y loop and the x loop is 1 because either there is
    a border area next to the halo -- then the first column or row is covered below in
    the border update -- or there is an outside border -- then the first column or row
    contains the boundary values. */

    for ( size_t y= 1; y < lh-1; y++ ) {

        const size_t x= 1;
        double* p_here= &gridlocalbegin[ y*lw+x ];
        double* p_east= &gridlocalbegin[ y*lw+x+1 ];
        double* p_west= &gridlocalbegin[ y*lw+x-1 ];
        double* p_north= &gridlocalbegin[ (y-1)*lw+x ];
        double* p_south= &gridlocalbegin[ (y+1)*lw+x ];

        for ( size_t x= 1; x < lw-1; x++ ) {

            double dtheta = ( *p_east + *p_west + *p_north + *p_south ) / 4.0 - *p_here ;
            *p_here += c * dtheta;
            res= std::max( res, std::fabs( dtheta ) );
            p_here++;
            p_east++;
            p_west++;
            p_north++;
            p_south++;
        }
    }
    MiniMonT::MiniMonRecord( 1, "smooth_inner", param );

    MiniMonT::MiniMonRecord( 0, "wait_outer ", param );

    // wait for async halo update
    level.halo.waitHalosAsync();

    MiniMonT::MiniMonRecord( 1, "wait_outer ", param );

    MiniMonT::MiniMonRecord( 0, "smooth_outer ", param );

    auto bend = level.halo.bend();
    // update border area
    for( auto it = level.halo.bbegin(); it != bend; ++it ) {

        double dtheta = ( it.valueAt(2) + it.valueAt(3) +
            it.valueAt(0) + it.valueAt(1) ) / 4.0 - *it;
        gridlocalbegin[ it.lpos() ] += c * dtheta;
        res= std::max( res, std::fabs( dtheta ) );
    }

    MiniMonT::MiniMonRecord( 1, "smooth_outer ", param );

    static dash::Array<double> residuals( dash::Team::All().size(), dash::BLOCKED );
    residuals.local[0]= res;

    MiniMonT::MiniMonRecord( 0, "sync_residual", param );
    residuals.barrier();

    if ( 0 != dash::myid() ) {

        dash::transform<double>(
            residuals.lbegin(), residuals.lend(), // first source
            residuals.begin(), // second source
            residuals.begin(), // destination
            dash::max<double>() );
    }
    residuals.barrier();
    MiniMonT::MiniMonRecord( 1, "sync_residual", param );

    MiniMonT::MiniMonRecord( 1, "smoothen", param );

    return residuals[0];
}


void v_cycle( vector<Level*>& levels, double epsilon= 0.01 ) {

    MiniMonT::MiniMonRecord( 0, "v_cycle" );

    writeToPNG( levels[0]->grid );

    for ( auto i= 1; i < levels.size(); i++ ) {

        double res= smoothen( *levels[i-1] );

        writeToPNG( levels[i-1]->grid );

        scaledown( *levels[i-1], *levels[i] );

        if ( 0 == dash::myid() ) {
            cout << "smoothen with residual " << res <<
                ", then scale down " << levels[i-1]->grid.extent(1) << "x" << levels[i-1]->grid.extent(1) <<
                " --> " << levels[i]->grid.extent(1)   << "x" << levels[i]->grid.extent(1) << endl;
        }

        writeToPNG( levels[i]->grid );
    }

    double residual= 1.0+epsilon;
    while ( residual > epsilon ) {

        residual= smoothen( *levels.back() );

        if ( 0 == dash::myid() ) {
            cout << "smoothen coarsest with residual " << residual << " > eps " << endl;
        }
        writeToPNG( levels.back()->grid );
    }

    for ( auto i= levels.size()-1; i > 0; i-- ) {

        scaleup( *levels[i], *levels[i-1] );

        writeToPNG( levels[i-1]->grid );

        double res= smoothen( *levels[i-1] );

        if ( 0 == dash::myid() ) {
            cout << "scale up " << levels[i]->grid.extent(1)   << "x" << levels[i]->grid.extent(1) <<
                " --> " << levels[i-1]->grid.extent(1) << "x" << levels[i-1]->grid.extent(1) <<
                ", then smoothen with residual " << res << endl;
        }

        writeToPNG( levels[i-1]->grid );
    }

    MiniMonT::MiniMonRecord( 1, "v_cycle" );
}


void smoothen_final( vector<Level*>& levels, double epsilon= 0.01 ) {

    MiniMonT::MiniMonRecord( 0, "smoothen_final" );

    double residual= 1.0+epsilon;
    while ( residual > epsilon ) {

        residual= smoothen( *levels.front() );

        if ( 0 == dash::myid() ) {
            cout << "smoothen finest with residual " << residual << " > eps " << endl;
        }
        writeToPNG( levels.front()->grid );
    }

    MiniMonT::MiniMonRecord( 1, "smoothen_final" );
}


int main( int argc, char* argv[] ) {

    dash::init(&argc, &argv);

    MiniMonT::MiniMonInit();

    MiniMonT::MiniMonRecord( 0, "main" );

    TeamSpecT teamspec( dash::Team::All().size(), 1 );
    teamspec.balance_extents();

    /* determine factors for width and height such that every unit has a power of two
    extent in every dimension and that the area is close to a square
    with aspect ratio \in [0,75,1,5] */

    uint32_t factor_y= teamspec.num_units(0);
    uint32_t factor_x= teamspec.num_units(1);

    uint32_t factor_max= std::max( factor_y, factor_x );
    while ( factor_y < 0.75 * factor_max ) { factor_y *= 2; }
    while ( factor_x < 0.75 * factor_max ) { factor_x *= 2; }

    constexpr uint32_t howmanylevels= 9;
    vector<Level*> levels;
    levels.reserve( howmanylevels );

    /* create all grid levels, starting with the finest and ending with 2x2 */
    for ( auto l= 0; l < howmanylevels-0; l++ ) {

        if ( 0 == dash::myid() ) {
            cout << "level " << l << " is " <<
                (1<<(howmanylevels-l))*factor_y << "x" << (1<<(howmanylevels-l))*factor_x << endl;
        }

        levels.push_back(
            new Level( (1<<(howmanylevels-l))*factor_y, (1<<(howmanylevels-l))*factor_x ,
            dash::Team::All(), teamspec ) );

        dash::barrier();
        /*
        cout << "unit " << dash::myid() << " : " <<
            levels.back()->grid.local.extent(1) << " x " <<
            levels.back()->grid.local.extent(0) << endl;
        dash::barrier();
        */

        if ( 0 == l ) {
            sanitycheck( levels[0]->grid );
        }

        /* disabled for now */
        initboundary( *levels[l] );
    }

    dash::barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here
    but we do it for demonstration in the PNG images */
    initgrid( levels[0]->grid );
    markunits( levels[0]->grid );

    writeToPNG( levels[0]->grid );

    markboundary( *levels[0] );

    writeToPNG( levels[0]->grid );

    dash::barrier();

    v_cycle( levels, 0.01 );
    smoothen_final( levels, 0.1 );

    MiniMonT::MiniMonRecord( 1, "main" );
    MiniMonT::MiniMonPrint( dash::myid(), dash::Team::All().size() );

    dash::finalize();
    return 0;
}
