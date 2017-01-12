#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cstdio>

#include <libdash.h>
#include <dash/experimental/HaloMatrix.h>

#define WITHPNGOUTPUT 1
#ifdef WITHPNGOUTPUT
#include <png++/png.hpp>

uint32_t filenumber= 0;
char filename[101];

#define WRITETOPNG(GRID) \
    dash::barrier(); \
    if ( 0 == dash::myid() ) { \
        snprintf( filename, 100, "image%04i.png", filenumber++ ); \
        topng( GRID , filename, 0.0, 10.0 ); \
    } \
    dash::barrier();

#else /* WITHPNGOUTPUT */

#define WRITETOPNG(GRID)
    
#endif /* WITHPNGOUTPUT */

using std::cout;
using std::endl;
using std::setw;
using std::vector;


dash::HaloSpec<2> halospec({{ {-1,1}, {-1,1} }});


class Level {

public:

    dash::Matrix<double,2> grid;
    dash::HaloMatrix< dash::Matrix<double,2>, dash::HaloSpec<2> > halo;
    
    Level( size_t w, size_t h, dash::Team& team, dash::TeamSpec<2> teamspec ) :
        grid( dash::SizeSpec<2>( h+2, w+2 ),
              dash::DistributionSpec<2>( dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
        halo( grid, halospec ) {}
    ~Level() {}
};


void sanitycheck( const dash::Matrix<double,2>& grid  ) {
    
    /* check if the sum of the local extents of the matrix blocks sum up to 
     * the global extents, abort otherwise */
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


void initgrid( dash::Matrix<double,2>& grid ) {

    if ( 0 != dash::myid() ) return;

    /* do initialization only from unit 0 with global accesses. 
    This is not the fastest way and a parallel routine would not be too 
    complicated but it should demonstrate the possibility to trade 
    convenience against performance. */

    size_t w= grid.extent(1);
    size_t h= grid.extent(0);

    cout << "assign values to finest level" << endl;

    for( size_t i = 1; i < h-1; ++i ) {
    for( size_t k = 1; k < w-1; ++k ) {
        grid[i][k] = 7.0;
    }
    }

    for( size_t i = h/5; i < h*3/4; ++i ) {
    for( size_t k = w*1/5; k < w*1/3; ++k ) {
        grid[i][k] = 3.0;
    }
    }
}


void setboundary( dash::Matrix<double,2>& grid ) {
    
    if ( 0 != dash::myid() ) return;
    
    /* do initialization only from unit 0 with global accesses. 
    This is not the fastest way and a parallel routine would not be too 
    complicated but it should demonstrate the possibility to trade 
    convenience against performance. */
    
    size_t w= grid.extent(1);
    size_t h= grid.extent(0);

    /* set entire boundary to a low value */
    for ( size_t x= 0; x < w; x++ ) {
        grid[0  ][x]= 5.0;
        grid[h-1][x]= 5.0;
    }
    for ( size_t y= 0; y < h; y++ ) {
        grid[y][0  ]= 5.0;
        grid[y][w-1]= 5.0;
    }

    /* set parts of boundary to a high value */
    for ( size_t x= 0; x < w/3; x++ ) {
        grid[0  ][x]= 10.0;
    }
    for ( size_t y= 0; y < h*3/5; y++ ) {
        grid[y][0  ]= 10.0;
    }

    /* set other parts of boundary to a very low value */
    for ( size_t y= h*1/5; y < h*2/5; y++ ) {
        grid[y][w-1]= 0.0;
    }
    for ( size_t y= h*3/5; y < h*4/5; y++ ) {
        grid[y][w-1]= 0.0;
    }
}


void scaledownboundary( const dash::Matrix<double,2>& fine, dash::Matrix<double,2>& coarse ) {

    assert( (coarse.extent(1)-2)*2 +2 == fine.extent(1) );
    assert( (coarse.extent(0)-2)*2 +2 == fine.extent(0) );

    size_t wc= coarse.local.extent(1);
    size_t hc= coarse.local.extent(0);
    size_t wf= fine.local.extent(1);
    size_t hf= fine.local.extent(0);
    std::array< long int, 2 > cornerc= coarse.pattern().global( {0,0} );

    if ( 0 == cornerc[0] ) {

        /* if boundary in front */
        size_t startx= ( 0 == cornerc[1] ) ? 1 : 0;
        for ( size_t x= 0; x < wc-1; x++ ) {
            coarse.local[0][startx+x]= 0.5 * ( fine.local[0][startx+2*x+0] + fine.local[0][startx+2*x+1] );
        }

    } else {
        
        /* else boundary in the end -- this is only true for the restructed case with 4 units */
        size_t startx= ( 0 == cornerc[1] ) ? 1 : 0;
        for ( size_t x= 0; x < wc-1; x++ ) {
            coarse.local[hc-1][startx+x]= 0.5 * ( fine.local[hf-1][startx+2*x+0] + fine.local[hf-1][startx+2*x+1] );
        }
    }

    if ( 0 == cornerc[1] ) {

        /* if boundary in front */
        size_t starty= ( 0 == cornerc[0] ) ? 1 : 0;
        for ( size_t y= 0; y < hc-1; y++ ) {
            coarse.local[starty+y][0]= 0.5 * ( fine.local[starty+2*y+0][0] + fine.local[starty+2*y+1][0] );
        }

    } else {
        
        /* else boundary in the end -- this is only true for the restructed case with 4 units */
        size_t starty= ( 0 == cornerc[0] ) ? 1 : 0;
        for ( size_t y= 0; y < hc-1; y++ ) {
            coarse.local[starty+y][wc-1]= 0.5 * ( fine.local[starty+2*y+0][wf-1] + fine.local[starty+2*y+1][wf-1] );
        }
    }
}


void scaledown( const dash::Matrix<double,2>& fine, dash::Matrix<double,2>& coarse ) {

    assert( (coarse.extent(1)-2)*2 +2 == fine.extent(1) );
    assert( (coarse.extent(0)-2)*2 +2 == fine.extent(0) );

    size_t w= coarse.local.extent(1);
    size_t h= coarse.local.extent(0);
    std::array< long int, 2 > cornerc= coarse.pattern().global( {0,0} );

    size_t startx= ( 0 == cornerc[1] ) ? 1 : 0;
    size_t starty= ( 0 == cornerc[0] ) ? 1 : 0;
    for ( size_t y= 0; y < h-1 ; y++ ) {
        for ( size_t x= 0; x < w-1 ; x++ ) {

            coarse.local[starty+y][startx+x]= 0.25 * (
                fine.local[starty+2*y  ][startx+2*x  ] +
                fine.local[starty+2*y  ][startx+2*x+1] +
                fine.local[starty+2*y+1][startx+2*x  ] +
                fine.local[starty+2*y+1][startx+2*x+1] );
        }
    }

    dash::barrier();
}


void scaleup( const dash::Matrix<double,2>& coarse, dash::Matrix<double,2>& fine ) {

    assert( (coarse.extent(1)-2)*2 +2 == fine.extent(1) );
    assert( (coarse.extent(0)-2)*2 +2 == fine.extent(0) );

    size_t w= coarse.local.extent(1);
    size_t h= coarse.local.extent(0);
    std::array< long int, 2 > cornerc= coarse.pattern().global( {0,0} );
    size_t startx= ( 0 == cornerc[1] ) ? 1 : 0;
    size_t starty= ( 0 == cornerc[0] ) ? 1 : 0;

    for ( size_t y= 0; y < h-1 ; y++ ) {
        for ( size_t x= 0; x < w-1 ; x++ ) {

            double t= coarse.local[starty+y][startx+x];
            fine.local[starty+2*y  ][startx+2*x  ]= t;
            fine.local[starty+2*y  ][startx+2*x+1]= t;
            fine.local[starty+2*y+1][startx+2*x  ]= t;
            fine.local[starty+2*y+1][startx+2*x+1]= t;
        }
    }

    dash::barrier();
}


double smoothen( dash::Matrix<double,2>& grid,
               dash::HaloMatrix< dash::Matrix<double,2>, dash::HaloSpec<2> >& halo ) {

    double res= 0.0;
    
    /// relaxation coeff.
    const double c= 1.0;

    /// begin pointer of local block, needed because halo border iterator is read-only
    auto gridlocalbegin= grid.lbegin();

    // async halo update
    halo.updateHalosAsync();

    // update inner

    size_t lw= grid.local.extent(1);
    size_t lh= grid.local.extent(0);

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

    // wait for async halo update
    halo.waitHalosAsync();

    // update border area
    for( auto it = halo.bbegin(); it != halo.bend(); ++it ) {

        auto core = *it;
        double dtheta = ( it.halo_value(1,-1) + it.halo_value(1,1) +
                          it.halo_value(0,-1) + it.halo_value(0,1) ) / 4.0 - *it;
        gridlocalbegin[ it.lpos() ] += c * dtheta;
        res= std::max( res, std::fabs( dtheta ) );
    }

    dash::Array<double> residuals( dash::Team::All().size(), dash::BLOCKED );
    residuals.local[0]= res;

    residuals.barrier();

    if ( 0 != dash::myid() ) {

        dash::transform<double>(
            residuals.lbegin(), residuals.lend(), // first source
            residuals.begin(), // second source
            residuals.begin(), // destination
            dash::max<double>() );
    }

    residuals.barrier();
    
    /* only unit 0 returns global residual */
    return residuals.local[0];
}


#ifdef WITHPNGOUTPUT
void topng( dash::Matrix<double,2>& grid, const char* filename, 
        const double minvalue= 0.0, const double maxvalue= 100.0 ) {

    if ( 0 != dash::myid() ) return;

    uint32_t minsize= 512;
    uint32_t f= 1;
    size_t w= grid.extent(1);
    size_t h= grid.extent(0);

    if ( w < minsize || h < minsize ) {
        
        f= ( w <= h ) ? w : h;
        f= 1+minsize/f;
    }
    
    png::image< png::rgb_pixel > image( h*f, w*f );
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

    image.write( filename );
}
#endif /* WITHPNGOUTPUT */


void v_cycle( vector<Level*>& levels, uint32_t inneriterations= 1 ) {

    WRITETOPNG( levels[0]->grid )

    for ( uint32_t i= 1; i < levels.size(); i++ ) {

        double res= smoothen( levels[i-1]->grid, levels[i-1]->halo );

        WRITETOPNG( levels[i-1]->grid )

        scaledown( levels[i-1]->grid, levels[i]->grid );

        if ( 0 == dash::myid() ) {
            cout << "smoothen with residual " << res << 
                ", then scale down " << levels[i-1]->grid.extent(1) << "x" << levels[i-1]->grid.extent(1) << 
                " --> " << levels[i]->grid.extent(1)   << "x" << levels[i]->grid.extent(1) << endl;
        }

        WRITETOPNG( levels[i]->grid )
    }

    /* do a fixed number of smoothing steps until the residual on the coarsest grid is 
    what you want on the finest level in the end. */
    for ( uint32_t i= 0; i < inneriterations; i++ ) {
        
        double res= smoothen( levels.back()->grid, levels.back()->halo );

        if ( 0 == dash::myid() ) {
            cout << "smoothen coarsest with residual " << res << endl;
        }
        WRITETOPNG( levels.back()->grid )
    }
    
    for ( uint32_t i= levels.size()-1; i > 0; i-- ) {

        scaleup( levels[i]->grid, levels[i-1]->grid );

        WRITETOPNG( levels[i-1]->grid )

        double res= smoothen( levels[i-1]->grid, levels[i-1]->halo );

        if ( 0 == dash::myid() ) {
            cout << "scale up " << levels[i]->grid.extent(1)   << "x" << levels[i]->grid.extent(1) << 
                " --> " << levels[i-1]->grid.extent(1) << "x" << levels[i-1]->grid.extent(1) << 
                ", then smoothen with residual " << res << endl;
        }

        WRITETOPNG( levels[i-1]->grid )
    }
}


int main( int argc, char* argv[] ) {

    dash::init(&argc, &argv);
    
    dash::TeamSpec<2> teamspec( dash::Team::All().size(), 1 );
    teamspec.balance_extents();

    uint32_t howmanylevels= 11;
    vector<Level*> levels;
    levels.reserve( howmanylevels );

    /* create all grid levels, starting with the finest and ending with 4x4 */
    for ( uint32_t l= 0; l < howmanylevels-1; l++ ) {

        if ( 0 == dash::myid() ) {
            cout << "level " << l << " is " << (1<<(howmanylevels-l))+2 << "x" << (1<<(howmanylevels-l))+2 << endl;
        }
        levels.push_back( new Level( (1<<(howmanylevels-l)), (1<<(howmanylevels-l)), dash::Team::All(), teamspec ) );
        
        if ( 0 == l ) {
            sanitycheck( levels[0]->grid );
            setboundary( levels[0]->grid );
        } else {
            scaledownboundary( levels[l-1]->grid, levels[l]->grid );
        }
    }
   
    dash::barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here 
    but we do it for demonstration in the PNG images */
    initgrid( levels[0]->grid );

    dash::barrier();
    
    v_cycle( levels, 10 );

    double res= smoothen( levels.front()->grid, levels.front()->halo );
    if ( 0 == dash::myid() ) {
        cout << "final residual " << res << endl;
    }

    WRITETOPNG( levels.front()->grid )
        
    for ( auto ll : levels ) {
        delete( ll );
    }

    dash::finalize();
}
