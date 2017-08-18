#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cstdio>

#include <libdash.h>
#include <dash/experimental/HaloMatrixWrapper.h>


#define WITHCSVOUTPUT 1
#ifdef WITHCSVOUTPUT

uint32_t filenumber= 0;
const std::string filename("image.csv");

#endif /* WITHCSVOUTPUT */


using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::vector;

using TeamSpecT = dash::TeamSpec<3>;
using MatrixT = dash::NArray<double,3>;
using StencilT = dash::experimental::Stencil<3>;
using StencilSpecT = dash::experimental::StencilSpec<3,6>;
using CycleSpecT = dash::experimental::CycleSpec<3>;
using HaloMatrixWrapperT = dash::experimental::HaloMatrixWrapper<MatrixT,StencilSpecT>;


const StencilSpecT stencil_spec({
    StencilT(-1, 0, 0), StencilT( 1, 0, 0),
    StencilT( 0,-1, 0), StencilT( 0, 1, 0),
    StencilT( 0, 0,-1), StencilT( 0, 0, 1)});

const CycleSpecT cycle_spec(
    dash::experimental::Cycle::FIXED,
    dash::experimental::Cycle::FIXED,
    dash::experimental::Cycle::FIXED );

struct Level {

    MatrixT grid;
    HaloMatrixWrapperT halo;

    Level( size_t d, size_t h, size_t w, dash::Team& team, TeamSpecT teamspec ) :
        grid( dash::SizeSpec<3>( d, h, w ),
            dash::DistributionSpec<3>( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
        halo( grid, stencil_spec, cycle_spec ) {}
    Level() = delete;

    ~Level() {}
};


/* Write out the entire state as a CSV which can be loaded and vizualized
with for example paraview. For convenience for vozualization, make the output
constant size regardless of the grid level. Thus interpolate coarser grids or
reduce finer ones. */

size_t resolutionForCSVd= 0;
size_t resolutionForCSVh= 0;
size_t resolutionForCSVw= 0;


#ifdef WITHCSVOUTPUT
void toCSV( const MatrixT& grid ) {

    if ( 0 != dash::myid() ) return;

    size_t d= grid.extent(0);
    size_t h= grid.extent(1);
    size_t w= grid.extent(2);

    std::ofstream csvfile;
    csvfile.open( filename + '.' + std::to_string(filenumber++) );

    csvfile << "z coord,y coord,x coord,heat" << "\n";

    /* Write according to the fixed grid size. If the real grid is even finer
    then pick any point that matches the coordinates. If the real grid is
    coarser, then repeat one value from the coarser grid multiple times */

    size_t divd= ( resolutionForCSVd > d ) ? resolutionForCSVd/d : 1;
    size_t divh= ( resolutionForCSVh > h ) ? resolutionForCSVh/h : 1;
    size_t divw= ( resolutionForCSVw > w ) ? resolutionForCSVw/w : 1;

    size_t muld= ( resolutionForCSVd < d ) ? d/resolutionForCSVd : 1;
    size_t mulh= ( resolutionForCSVh < h ) ? h/resolutionForCSVh : 1;
    size_t mulw= ( resolutionForCSVw < w ) ? w/resolutionForCSVw : 1;

    for ( size_t z = 0; z < resolutionForCSVd; ++z ) {
        for ( size_t y = 0; y < resolutionForCSVh; ++y ) {
            for ( size_t x = 0; x < resolutionForCSVw; ++x ) {

                csvfile << x << "," << y << "," << z << "," <<
                    (double) grid[muld*z/divd][mulh*y/divh][mulw*x/divw] << "\n";
            }
        }
    }

    csvfile.close();
}
#endif /* WITHCSVOUTPUT */

void writeToCsv( const MatrixT& grid ) {
#ifdef WITHCSVOUTPUT
    grid.barrier();

    if ( 0 == dash::myid() )
        toCSV( grid );

    grid.barrier();
#endif /* WITHCSVOUTPUT */
}


void sanitycheck( const MatrixT& grid  ) {

    /* check if the sum of the local extents of the matrix blocks sum up to
    the global extents, abort otherwise */
    dash::Array<size_t> sums( dash::Team::All().size(), dash::BLOCKED );
    sums.local[0]= grid.local.extent(0) * grid.local.extent(1) * grid.local.extent(2);
    sums.barrier();

    if ( 0 != dash::myid() ) {

        dash::transform<size_t>(
            sums.lbegin(), sums.lend(), // first source
            sums.begin(), // second source
            sums.begin(), // destination
            dash::plus<size_t>() );
    }
    sums.barrier();

    if ( ( (size_t) sums[0] ) != grid.extent(0) * grid.extent(1) * grid.extent(2) ) {

        if ( 0 == dash::myid() ) {
            cout << "ERROR: size mismatch: global size is " <<
                grid.extent(0) * grid.extent(1) * grid.extent(2) << " == " <<
                grid.extent(0) << "x" << grid.extent(1) << "x" << grid.extent(2) <<
                " but sum of local sizes is " << (size_t) sums[0] << std::flush << endl;
        }

        dash::finalize();
        exit(0);
    }
}


void initgrid( MatrixT& grid ) {

    /* init to 0 first! */
    dash::fill( grid.begin(), grid.end(), 5.0 );

    /* do initialization only from unit 0 with global accesses.
    This is not the fastest way and a parallel routine would not be too
    complicated but it should demonstrate the possibility to trade
    convenience against performance. */

    if ( 0 == dash::myid() ) {


        size_t w= grid.extent(2);
        size_t h= grid.extent(1);
        size_t d= grid.extent(0);

        for ( size_t i = d*1/5; i < d*4/5; ++i ) {
            for ( size_t j = h/4; j < h/2; ++j ) {
                for ( size_t k = w*4/6; k < w*5/6; ++k ) {
                    grid[i][j][k] = 8.0;
                }
            }
        }

        for ( size_t i = d*1/5; i < d*4/5; ++i ) {
            for ( size_t j = h/5; j < h*3/4; ++j ) {
                for ( size_t k = w*1/5; k < w*1/3; ++k ) {
                    grid[i][j][k] = 2.0;
                }
            }
        }

    }

    grid.barrier();
}


void initboundary( Level& level ) {

    double gw= level.grid.extent(2);
    double gh= level.grid.extent(1);
    double gd= level.grid.extent(0);

    /* with this way of setting the boundary conditions on every level separatel (in
    contrast to scaling it down step by step) one needs to make sure that the boundary
    values are continuous (kind of). Otherwise the coarsened version of the boundary
    values will have jumps in different places (instead of jumps being smoothed out
    by coarsening the boundary).
    Well, we should really think about another way to init the boundary values! */

    level.halo.setFixedHalos( [gd,gh,gw]( const std::array<dash::default_index_t,3>& coords ) {

        dash::default_index_t z= coords[0];
        dash::default_index_t y= coords[1];
        dash::default_index_t x= coords[2];
        double ret= 0;

        /* for simplicity make every side uniform */

        if ( -1 == z ) {

            ret= 10.0;

        } else if ( gd == z ) {

            ret= 0.0;

        } else if ( -1 == y ) {

            ret= 9.0;

        } else if ( gh == y ) {

            ret= 1.0;

        } else if ( -1 == x ) {

            ret= 8.0;

        } else if ( gw == x ) {

            ret= 2.0;
        }

        return ret;
    });

}


void markunits( MatrixT& grid ) {

    /* Mark unit bordery by setting the first local rows and columns */

    size_t w= grid.local.extent(2);
    size_t h= grid.local.extent(1);
    size_t d= grid.local.extent(0);

    for ( size_t i = 0; i < d; ++i ) {
        for ( size_t j = 0; j < h; ++j ) {
            grid.local[i][j][0] = 9.0;
        }
    }

    for ( size_t i = 0; i < d; ++i ) {
        for ( size_t k = 0; k < w; ++k ) {
            grid.local[i][0][k] = 9.0;
        }
    }

    for ( size_t j = 0; j < h; ++j ) {
        for ( size_t k = 0; k < w; ++k ) {
            grid.local[0][j][k] = 9.0;
        }
    }

}


/* currently not needed */
void scaledownboundary( const MatrixT& fine, MatrixT& coarse ) {

    assert( (coarse.extent(2)-2)*2 +2 == fine.extent(2) );
    assert( (coarse.extent(1)-2)*2 +2 == fine.extent(1) );
    assert( (coarse.extent(0)-2)*2 +2 == fine.extent(0) );

    size_t wc= coarse.local.extent(2);
    size_t hc= coarse.local.extent(1);
    size_t dc= coarse.local.extent(0);
    size_t wf= fine.local.extent(2);
    size_t hf= fine.local.extent(1);
    size_t df= fine.local.extent(0);
    std::array< long int, 3 > cornerc= coarse.pattern().global( {0,0,0} );

#if 0
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
#endif
}


void scaledown( Level& fine, Level& coarse ) {

    assert( coarse.grid.extent(2) * 2 == fine.grid.extent(2) );
    assert( coarse.grid.extent(1) * 2 == fine.grid.extent(1) );
    assert( coarse.grid.extent(0) * 2 == fine.grid.extent(0) );

    const std::array< long unsigned int, 3 > extentc= coarse.grid.pattern().local_extents();
    const std::array< long signed int, 3 >   cornerc= coarse.grid.pattern().global( {0,0,0} );
    const std::array< long unsigned int, 3 > extentf= fine.grid.pattern().local_extents();
    const std::array< long signed int, 3 >   cornerf= fine.grid.pattern().global( {0,0,0} );

    assert( cornerc[0] * 2 == cornerf[0] );
    assert( cornerc[1] * 2 == cornerf[1] );
    assert( cornerc[2] * 2 == cornerf[2] );

    assert( 0 == cornerc[0] %2 );
    assert( 0 == cornerc[1] %2 );
    assert( 0 == cornerc[2] %2 );

    assert( extentc[0] * 2 == extentf[0] );
    assert( extentc[1] * 2 == extentf[1] );
    assert( extentc[2] * 2 == extentf[2] );

    for ( size_t z= 0; z < extentc[0]; z++ ) {
        for ( size_t y= 0; y < extentc[1]; y++ ) {
            for ( size_t x= 0; x < extentc[2]; x++ ) {

                coarse.grid.local[z][y][x]= 1.0 / 8.0 * (
                    fine.grid.local[2*z  ][2*y  ][2*x  ] +
                    fine.grid.local[2*z  ][2*y  ][2*x+1] +
                    fine.grid.local[2*z  ][2*y+1][2*x  ] +
                    fine.grid.local[2*z  ][2*y+1][2*x+1] +
                    fine.grid.local[2*z+1][2*y  ][2*x  ] +
                    fine.grid.local[2*z+1][2*y  ][2*x+1] +
                    fine.grid.local[2*z+1][2*y+1][2*x  ] +
                    fine.grid.local[2*z+1][2*y+1][2*x+1] );
            }
        }
    }

}


void scaleup( Level& coarse, Level& fine ) {

    assert( coarse.grid.extent(2) * 2 == fine.grid.extent(2) );
    assert( coarse.grid.extent(1) * 2 == fine.grid.extent(1) );
    assert( coarse.grid.extent(0) * 2 == fine.grid.extent(0) );

    const std::array< long unsigned int, 3 > extentc= coarse.grid.pattern().local_extents();
    const std::array< long signed int, 3 >   cornerc= coarse.grid.pattern().global( {0,0,0} );
    const std::array< long unsigned int, 3 > extentf= fine.grid.pattern().local_extents();
    const std::array< long signed int, 3 >   cornerf= fine.grid.pattern().global( {0,0,0} );

    assert( cornerc[0] * 2 == cornerf[0] );
    assert( cornerc[1] * 2 == cornerf[1] );
    assert( cornerc[2] * 2 == cornerf[2] );

    assert( 0 == cornerc[0] %2 );
    assert( 0 == cornerc[1] %2 );
    assert( 0 == cornerc[2] %2 );

    assert( extentc[0] * 2 == extentf[0] );
    assert( extentc[1] * 2 == extentf[1] );
    assert( extentc[2] * 2 == extentf[2] );

    for ( size_t z= 0; z < extentc[0] ; z++ ) {
        for ( size_t y= 0; y < extentc[1] ; y++ ) {
            for ( size_t x= 0; x < extentc[2]; x++ ) {

                double t= coarse.grid.local[z][y][x];
                fine.grid.local[2*z  ][2*y  ][2*x  ]= t;
                fine.grid.local[2*z  ][2*y  ][2*x+1]= t;
                fine.grid.local[2*z  ][2*y+1][2*x  ]= t;
                fine.grid.local[2*z  ][2*y+1][2*x+1]= t;
                fine.grid.local[2*z+1][2*y  ][2*x  ]= t;
                fine.grid.local[2*z+1][2*y  ][2*x+1]= t;
                fine.grid.local[2*z+1][2*y+1][2*x  ]= t;
                fine.grid.local[2*z+1][2*y+1][2*x+1]= t;
            }
        }
    }
}


double smoothen( Level& level ) {

    double res= 0.0;

    /// relaxation coeff.
    const double c= 1.0;

    /* this barrier is there so that iterations are synchronized across all
    units. Otherwise some overtak others esp. on the very small grids. */
    level.grid.barrier();

    /// begin pointer of local block, needed because halo border iterator is read-only
    auto gridlocalbegin= level.grid.lbegin();

    // async halo update
    level.halo.updateHalosAsync();

    // update inner

    size_t lw= level.grid.local.extent(2);
    size_t lh= level.grid.local.extent(1);
    size_t ld= level.grid.local.extent(0);

    /* the start value for both, the y loop and the x loop is 1 because either there is
    a border area next to the halo -- then the first column or row is covered below in
    the border update -- or there is an outside border -- then the first column or row
    contains the boundary values. */

    for ( size_t z= 1; z < ld-1; z++ ) {
        for ( size_t y= 1; y < lh-1; y++ ) {

            const size_t x= 1;
            double* p_here=  &gridlocalbegin[ ( (z)*lh+(y))*lw+(x) ];
            double* p_east=  &gridlocalbegin[ ( (z)*lh+(y))*lw+(x+1) ];
            double* p_west=  &gridlocalbegin[ ( (z)*lh+(y))*lw+(x-1) ];
            double* p_north= &gridlocalbegin[ ( (z)*lh+(y-1))*lw+(x) ];
            double* p_south= &gridlocalbegin[ ( (z)*lh+(y+1))*lw+(x) ];
            double* p_up=    &gridlocalbegin[ ( (z-1)*lh+(y))*lw+(x) ];
            double* p_down=  &gridlocalbegin[ ( (z+1)*lh+(y))*lw+(x) ];

            for ( size_t x= 1; x < lw-1; x++ ) {

                double dtheta = ( *p_east + *p_west + *p_north + *p_south + *p_up + *p_down ) / 6.0 - *p_here ;
                *p_here += c * dtheta;
                res= std::max( res, std::fabs( dtheta ) );
                p_here++;
                p_east++;
                p_west++;
                p_north++;
                p_south++;
                p_up++;
                p_down++;
            }
        }
    }

    // wait for async halo update
    level.halo.waitHalosAsync();

    auto bend = level.halo.bend();
    // update border area
    for( auto it = level.halo.bbegin(); it != bend; ++it ) {

        double dtheta = ( it.valueAt(0) + it.valueAt(1) +
            it.valueAt(2) + it.valueAt(3) + it.valueAt(4) + it.valueAt(5) ) / 6.0 - *it;
        gridlocalbegin[ it.lpos() ] += c * dtheta;
        res= std::max( res, std::fabs( dtheta ) );
    }

    static dash::Array<double> residuals( dash::Team::All().size(), dash::BLOCKED );
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

    return residuals[0];
}


void v_cycle( vector<Level*>& levels, double epsilon= 0.01 ) {

    writeToCsv( levels[0]->grid );

    for ( auto i= 1; i < levels.size(); i++ ) {

        double res= smoothen( *levels[i-1] );

        writeToCsv( levels[i-1]->grid );

        scaledown( *levels[i-1], *levels[i] );

        if ( 0 == dash::myid() ) {
            cout << "smoothen with residual " << res <<
                ", then scale down " << levels[i-1]->grid.extent(2) << "x" <<
                levels[i-1]->grid.extent(1) << "x"<< levels[i-1]->grid.extent(0) <<
                " --> " << levels[i]->grid.extent(2) << "x" <<
                levels[i]->grid.extent(1) << "x" << levels[i]->grid.extent(0) << endl;
        }

        writeToCsv( levels[i]->grid );
    }

    double residual= 1.0+epsilon;
    while ( residual > epsilon ) {

        residual= smoothen( *levels.back() );

        if ( 0 == dash::myid() ) {
            cout << "smoothen coarsest with residual " << residual << endl;
        }
        writeToCsv( levels.back()->grid );
    }

    for ( auto i= levels.size()-1; i > 0; i-- ) {

        scaleup( *levels[i], *levels[i-1] );

        writeToCsv( levels[i-1]->grid );

        double res= smoothen( *levels[i-1] );

        if ( 0 == dash::myid() ) {
            cout << "scale up " << levels[i]->grid.extent(2) << "x" <<
                levels[i]->grid.extent(1) << "x" << levels[i]->grid.extent(0) <<
                " --> " << levels[i-1]->grid.extent(2) << "x" <<
                levels[i-1]->grid.extent(1) << "x" << levels[i-1]->grid.extent(0) <<
                ", then smoothen with residual " << res << endl;
        }

        writeToCsv( levels[i-1]->grid );
    }
}


void smoothen_final( vector<Level*>& levels, double epsilon= 0.01 ) {

    double residual= 1.0+epsilon;
    while ( residual > epsilon ) {

        residual= smoothen( *levels.front() );

        if ( 0 == dash::myid() ) {
            cout << "smoothen finest with residual " << residual << endl;
        }
        writeToCsv( levels.front()->grid );
    }
}


int main( int argc, char* argv[] ) {

    dash::init(&argc, &argv);

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    /* determine factors for width and height such that every unit has a power of two
    extent in every dimension and that the area is close to a square
    with aspect ratio \in [0,75,1,5] */

    uint32_t factor_z= teamspec.num_units(0);
    uint32_t factor_y= teamspec.num_units(1);
    uint32_t factor_x= teamspec.num_units(2);

    uint32_t factor_max= factor_z;
    factor_max= std::max( factor_max, factor_y );
    factor_max= std::max( factor_max, factor_x );
    while ( factor_z < 0.75 * factor_max ) { factor_z *= 2; }
    while ( factor_y < 0.75 * factor_max ) { factor_y *= 2; }
    while ( factor_x < 0.75 * factor_max ) { factor_x *= 2; }

    constexpr uint32_t howmanylevels= 5;
    vector<Level*> levels;
    levels.reserve( howmanylevels );

    resolutionForCSVd= ( 1<<(howmanylevels) ) * factor_z;
    resolutionForCSVh= ( 1<<(howmanylevels) ) * factor_y;
    resolutionForCSVw= ( 1<<(howmanylevels) ) * factor_x;

    /* create all grid levels, starting with the finest and ending with 2x2 */
    for ( auto l= 0; l < howmanylevels-0; l++ ) {

        if ( 0 == dash::myid() ) {
            cout << "level " << l << " is " <<
                (1<<(howmanylevels-l))*factor_z << "x" <<
                (1<<(howmanylevels-l))*factor_y << "x" <<
                (1<<(howmanylevels-l))*factor_x << endl;
        }

        levels.push_back(
            new Level( (1<<(howmanylevels-l))*factor_z ,
                        (1<<(howmanylevels-l))*factor_y ,
                        (1<<(howmanylevels-l))*factor_x ,
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

    writeToCsv( levels[0]->grid );

    writeToCsv( levels[0]->grid );

    dash::barrier();

    v_cycle( levels, 0.01 );
    smoothen_final( levels, 0.1 );

    dash::finalize();
    return 0;
}
