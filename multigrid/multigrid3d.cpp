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

#define WITHCSVOUTPUT 1
#ifdef WITHCSVOUTPUT

uint32_t filenumber= 0;

#endif /* WITHCSVOUTPUT */


/* for MiniMonT */
std::vector<MiniMonT>* MiniMonT::tape;
std::chrono::time_point<std::chrono::high_resolution_clock> MiniMonT::inittime;


using std::cout;
using std::setfill;
using std::setw;
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
with for example paraview. For convenience for vizualization, make the output
constant size regardless of the grid level. Thus interpolate coarser grids or
reduce finer ones. */

size_t resolutionForCSVd= 0;
size_t resolutionForCSVh= 0;
size_t resolutionForCSVw= 0;


/* write out the current state of the given grid as CSV so that Paraview can
read and visualize it as a structured grid. Use a fixed size for the grid as
defined by the previous three global variables. So an animations of the
multigrid procedure is possible. */
void writeToCsv( const MatrixT& grid ) {

#ifdef WITHCSVOUTPUT

    grid.barrier();

    std::array< long int, 3 > corner= grid.pattern().global( {0,0,0} );

    size_t d= grid.extent(0);
    size_t h= grid.extent(1);
    size_t w= grid.extent(2);

    size_t dl= grid.local.extent(0);
    size_t hl= grid.local.extent(1);
    size_t wl= grid.local.extent(2);

    std::ofstream csvfile;
    csvfile.open( "image_unit" + std::to_string(dash::myid()) +
        ".csv." + std::to_string(filenumber++) );

    if ( 0 == dash::myid() ) {
        csvfile << " z coord, y coord, x coord, heat" << "\n";
    }

    /* Write according to the fixed grid size. If the real grid is even finer
    then pick any point that matches the coordinates. If the real grid is
    coarser, then repeat one value from the coarser grid multiple times */

    size_t divd= ( resolutionForCSVd > d ) ? resolutionForCSVd/d : 1;
    size_t divh= ( resolutionForCSVh > h ) ? resolutionForCSVh/h : 1;
    size_t divw= ( resolutionForCSVw > w ) ? resolutionForCSVw/w : 1;
    size_t muld= ( resolutionForCSVd < d ) ? d/resolutionForCSVd : 1;
    size_t mulh= ( resolutionForCSVh < h ) ? h/resolutionForCSVh : 1;
    size_t mulw= ( resolutionForCSVw < w ) ? w/resolutionForCSVw : 1;

    /* this should divide witout remainder */
    size_t localResolutionForCSVd= resolutionForCSVd * dl / d;
    size_t localResolutionForCSVh= resolutionForCSVh * hl / h;
    size_t localResolutionForCSVw= resolutionForCSVw * wl / w;

    corner[0]= corner[0] * resolutionForCSVd / d;
    corner[1]= corner[1] * resolutionForCSVh / h;
    corner[2]= corner[2] * resolutionForCSVw / w;

    for ( size_t z = 0 ; z < localResolutionForCSVd; ++z ) {
        for ( size_t y = 0; y < localResolutionForCSVh; ++y ) {
            for ( size_t x = 0; x < localResolutionForCSVw; ++x ) {

                csvfile << setfill('0') << setw(4) << corner[0]+z << "," <<
                    setfill('0') << setw(4) << corner[1]+y << "," <<
                    setfill('0') << setw(4) << corner[2]+x << "," <<
                    (double) grid.local[ muld*z/divd ][ mulh*y/divh ][ mulw*x/divw ] << "\n";
            }
        }
    }

    csvfile.close();
    grid.barrier();
#endif /* WITHCSVOUTPUT */
}


/* write out the grid in its actual size */
void writeToCsvFullGrid( const MatrixT& grid ) {

#ifdef WITHCSVOUTPUT

    grid.barrier();

    std::array< long int, 3 > corner= grid.pattern().global( {0,0,0} );

    size_t d= grid.extent(0);
    size_t h= grid.extent(1);
    size_t w= grid.extent(2);

    size_t dl= grid.local.extent(0);
    size_t hl= grid.local.extent(1);
    size_t wl= grid.local.extent(2);

    std::ofstream csvfile;
    csvfile.open( "image_unit" + std::to_string(dash::myid()) +
        ".csv." + std::to_string(filenumber++) );

    if ( 0 == dash::myid() ) {
        csvfile << " z coord, y coord, x coord, heat" << "\n";
    }

    for ( size_t z = 0 ; z < dl; ++z ) {
        for ( size_t y = 0; y < hl; ++y ) {
            for ( size_t x = 0; x < wl; ++x ) {

                csvfile << setfill('0') << setw(4) << corner[0]+z << "," <<
                    setfill('0') << setw(4) << corner[1]+y << "," <<
                    setfill('0') << setw(4) << corner[2]+x << "," <<
                    (double) grid.local[ z ][ y ][ x ] << "\n";
            }
        }
    }

    csvfile.close();
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

    /* not strictly necessary but it also avoids NAN values */
    dash::fill( grid.begin(), grid.end(), 5.0 );

    grid.barrier();
}

/* simple boundary settings where 3 sides are hot and 3 are cold */
void initboundary_3hot3cold( Level& level ) {

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


/* alternative boundary settings, the top and bottom planes have a hot circle
in the middle */
// void initboundary_2circles( Level& level ) {
void initboundary( Level& level ) {

    double gw= level.grid.extent(2);
    double gh= level.grid.extent(1);
    double gd= level.grid.extent(0);

    /* This way of setting boundaries uses subsampling on the top and bottom
    planes to determine the border values. This is another logical way that
    may be convenient sometimes. It guarantees that the boundary values on all
    the levels match.
    All other sides are constant at 0.0 degrees. The top an bottom circles are
    hot with 10.0 degrees. */

    level.halo.setFixedHalos( [gd,gh,gw]( const std::array<dash::default_index_t,3>& coords ) {

        dash::default_index_t z= coords[0];
        dash::default_index_t y= coords[1];
        dash::default_index_t x= coords[2];

        double ret= 0.0;

        /* for simplicity make every side uniform */

        if ( -1 == z || gd == z ) {

            /* radius differs on top and bottom plane */
            //double r= ( -1 == z ) ? 0.4 : 0.3;
            double r= 0.4;
            double r2= r*r;

            double lowvalue= 0.0;
            double highvalue= 10.0;

            double midx= 0.5;
            double midy= 0.5;

            /* At entry (x/gw,y/gh) we sample the
            rectangle [ x/gw,(x+1)/gw ) x [ y/gw, (y+1)/gh ) with m² points. */
            int32_t m= 5;
            int32_t m2= m*m;

            double sum= 0.0;
            double weight= 0.0;

            for ( int32_t iy= -m+1; iy < m; iy++ ) {
                for ( int32_t ix= -m+1; ix < m; ix++ ) {

                    double sx= (x+ix/m)/(gw-1);
                    double sy= (y+iy/m)/(gh-1);

                    double d2= (sx-midx)*(sx-midx) + (sy-midy)*(sy-midy);
                    sum += ( d2 <= r2 ) ? highvalue : lowvalue;
                    weight += 1.0;
                }
            }
            ret = sum / weight;
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

#if 0

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
#endif


void scaledown( Level& fine, Level& coarse ) {

    uint64_t param= fine.grid.local.extent(0)*fine.grid.local.extent(1)*fine.grid.local.extent(2);
    MiniMonT::MiniMonRecord( 0, "scaledown", param );

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

            const size_t x= 0;
            double* p_coarse= &coarse.grid.local[z][y][x];

            double* p_000= &fine.grid.local[2*z+0][2*y+0][2*x];
            double* p_010= &fine.grid.local[2*z+0][2*y+1][2*x];
            double* p_100= &fine.grid.local[2*z+1][2*y+0][2*x];
            double* p_110= &fine.grid.local[2*z+1][2*y+1][2*x];

            for ( size_t x= 0; x < extentc[2]; x++ ) {

                *p_coarse= 1.0 / 8.0 * (
                    *(p_000+0) + *(p_000+1) +
                    *(p_010+0) + *(p_010+1) +
                    *(p_100+0) + *(p_100+1) +
                    *(p_110+0) + *(p_110+1) );

                p_coarse++;
                p_000 += 2;
                p_010 += 2;
                p_100 += 2;
                p_110 += 2;
            }
        }
    }

    MiniMonT::MiniMonRecord( 1, "scaledown", param );
}


void scaleup( Level& coarse, Level& fine ) {

    uint64_t param= coarse.grid.local.extent(0)*coarse.grid.local.extent(1)*coarse.grid.local.extent(2);
    MiniMonT::MiniMonRecord( 0, "scaleup", param );

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

            const size_t x= 0;
            double* p_coarse= &coarse.grid.local[z][y][x];

            double* p_000= &fine.grid.local[2*z+0][2*y+0][2*x];
            double* p_010= &fine.grid.local[2*z+0][2*y+1][2*x];
            double* p_100= &fine.grid.local[2*z+1][2*y+0][2*x];
            double* p_110= &fine.grid.local[2*z+1][2*y+1][2*x];

            for ( size_t x= 0; x < extentc[2]; x++ ) {

                *(p_000+0)= *p_coarse;
                *(p_000+1)= *p_coarse;
                *(p_010+0)= *p_coarse;
                *(p_010+1)= *p_coarse;
                *(p_100+0)= *p_coarse;
                *(p_100+1)= *p_coarse;
                *(p_110+0)= *p_coarse;
                *(p_110+1)= *p_coarse;

                p_coarse++;
                p_000 += 2;
                p_010 += 2;
                p_100 += 2;
                p_110 += 2;
            }
        }
    }

    MiniMonT::MiniMonRecord( 1, "scaleup", param );
}


/* the parallel global residual is returned as a return parameter, but only
if it is not NULL because then the expensive parallel reduction is just avoided */
void smoothen( Level& level, uint32_t iter, double* residual_ret= NULL ) {

    uint64_t param= level.grid.local.extent(0)*level.grid.local.extent(1)*level.grid.local.extent(2);
    MiniMonT::MiniMonRecord( 0, "smoothen", param );

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

    MiniMonT::MiniMonRecord( 0, "smooth_inner", param );
    // update inner

    size_t lw= level.grid.local.extent(2);
    size_t lh= level.grid.local.extent(1);
    size_t ld= level.grid.local.extent(0);

    /* the start value for both, the y loop and the x loop is 1 because either there is
    a border area next to the halo -- then the first column or row is covered below in
    the border update -- or there is an outside border -- then the first column or row
    contains the boundary values. */

    if ( 0 == iter % 2 ) {

        for ( size_t z= 1; z < ld-1; z++ ) {
            for ( size_t y= 1; y < lh-1; y++ ) {

                /* this should eventually be done with Alpaka or Kokkos to look
                much nicer but still be fast */

                const size_t x= 1;
                double* p_here=  &level.grid.local[z  ][y  ][x  ];
                double* p_east=  &level.grid.local[z  ][y  ][x+1];
                double* p_west=  &level.grid.local[z  ][y  ][x-1];
                double* p_north= &level.grid.local[z  ][y+1][x  ];
                double* p_south= &level.grid.local[z  ][y-1][x  ];
                double* p_up=    &level.grid.local[z+1][y  ][x  ];
                double* p_down=  &level.grid.local[z-1][y  ][x  ];

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
    } else {
        for ( size_t z= ld-2; z >=1; z-- ) {
            for ( size_t y= lh-2; y >= 1; y-- ) {

                /* this should eventually be done with Alpaka or Kokkos to look
                much nicer but still be fast */

                const size_t x= lw-2;
                double* p_here=  &level.grid.local[z  ][y  ][x  ];
                double* p_east=  &level.grid.local[z  ][y  ][x+1];
                double* p_west=  &level.grid.local[z  ][y  ][x-1];
                double* p_north= &level.grid.local[z  ][y+1][x  ];
                double* p_south= &level.grid.local[z  ][y-1][x  ];
                double* p_up=    &level.grid.local[z+1][y  ][x  ];
                double* p_down=  &level.grid.local[z-1][y  ][x  ];

                for ( size_t x= lw-2; x >= 1; x-- ) {

                    double dtheta = ( *p_east + *p_west + *p_north + *p_south + *p_up + *p_down ) / 6.0 - *p_here ;
                    *p_here += c * dtheta;
                    res= std::max( res, std::fabs( dtheta ) );
                    p_here--;
                    p_east--;
                    p_west--;
                    p_north--;
                    p_south--;
                    p_up--;
                    p_down--;
                }
            }
        }
    }

    MiniMonT::MiniMonRecord( 1, "smooth_inner", param );

    MiniMonT::MiniMonRecord( 0, "smooth_wait", param );

    // wait for async halo update
    level.halo.waitHalosAsync();

    MiniMonT::MiniMonRecord( 1, "smooth_wait", param );

    MiniMonT::MiniMonRecord( 0, "smooth_outer", param );

    auto bend = level.halo.bend();
    // update border area
    for( auto it = level.halo.bbegin(); it != bend; ++it ) {

        double dtheta = ( it.valueAt(0) + it.valueAt(1) +
            it.valueAt(2) + it.valueAt(3) + it.valueAt(4) + it.valueAt(5) ) / 6.0 - *it;
        gridlocalbegin[ it.lpos() ] += c * dtheta;
        res= std::max( res, std::fabs( dtheta ) );
    }

    MiniMonT::MiniMonRecord( 1, "smooth_outer", param );

    if ( NULL != residual_ret ) {

        MiniMonT::MiniMonRecord( 0, "smooth_residuals", param );

        /*static*/ dash::Array<double> residuals( level.grid.team().size(), dash::BLOCKED, level.grid.team() );
        residuals.local[0]= res;

    /*
        residuals.barrier();

        if ( 0 != dash::myid() ) {

            dash::transform<double>(
                residuals.lbegin(), residuals.lend(), // first source
                residuals.begin(), // second source
                residuals.begin(), // destination
                dash::max<double>() );
        }
        residuals.barrier();
    */

        residuals.barrier();

        *residual_ret= 0.0;
        if ( 0 == dash::myid() ) {
            for ( const double& residual : residuals ) {
                *residual_ret= std::max( *residual_ret, residual );
            }
            residuals.local[0]= *residual_ret;
        }
        residuals.barrier();

        *residual_ret= residuals[0];

        MiniMonT::MiniMonRecord( 1, "smooth_residuals", param );
    }

    MiniMonT::MiniMonRecord( 1, "smoothen", param );
}


/* recursive version */
void v_cycle( vector<Level*>::const_iterator it, vector<Level*>::const_iterator itend,
        uint32_t numiter, double epsilon= 0.01 ) {

    if ( 0 == dash::myid() ) {
        cout << "v-cycle on  " <<
                    (*it)->grid.extent(2) << "x" <<
                    (*it)->grid.extent(1) << "x" <<
                    (*it)->grid.extent(0) << endl;
    }

    vector<Level*>::const_iterator itnext( it );
    itnext++;

    /* reached end of recursion? */
    if ( itend == itnext ) {

        /* smoothen completely  */

        double residual= 1.0+epsilon;
        uint32_t j= 0;
        while ( residual > epsilon ) {

            /* need global residual for iteration count */
            smoothen( **it, j++, &residual );

            if ( 0 == dash::myid() ) {
                cout << j << ": smoothen coarsest with residual " << residual << endl;
            }
        }
        writeToCsv( (*it)->grid );

        return;
    }

    /* normal recursion */

    /* smoothen somewhat with fixed number of iterations */
    for ( uint32_t j= 0; j < numiter; j++ ) {
        smoothen( **it, j, NULL );
    }

    writeToCsv( (*it)->grid );

    /* scale down */
    if ( 0 == dash::myid() ) {
        cout << "scale down " <<
            (*it)->grid.extent(2) << "x" <<
            (*it)->grid.extent(1) << "x" <<
            (*it)->grid.extent(0) <<
            " --> " <<
            (*itnext)->grid.extent(2) << "x" <<
            (*itnext)->grid.extent(1) << "x" <<
            (*itnext)->grid.extent(0) << endl;
    }
    scaledown( **it, **itnext );
    writeToCsv( (*it)->grid );

    /* recurse  */
    v_cycle( itnext, itend, numiter, epsilon );

    /* scale up */
    if ( 0 == dash::myid() ) {
        cout << "scale up " <<
            (*itnext)->grid.extent(2) << "x" <<
            (*itnext)->grid.extent(1) << "x" <<
            (*itnext)->grid.extent(0) <<
            " --> " <<
            (*it)->grid.extent(2) << "x" <<
            (*it)->grid.extent(1) << "x" <<
            (*it)->grid.extent(0) << endl;
    }
    scaleup( **itnext, **it );
    writeToCsv( (*it)->grid );

    /* smoothen somewhat with fixed number of iterations */
    for ( uint32_t j= 0; j < numiter; j++ ) {
        smoothen( **it, j, NULL );
    }
    writeToCsv( (*it)->grid );
}


void smoothen_final( vector<Level*>& levels, double epsilon= 0.01 ) {

    MiniMonT::MiniMonRecord( 0, "smoothfinal" );

    double residual= 1.0+epsilon;
    uint32_t j= 0;
    while ( residual > epsilon ) {

        smoothen( *levels.front(), j++, &residual );

        if ( 0 == dash::myid() ) {
            cout << j << ": smoothen finest with residual " << residual << endl;
        }
        //writeToCsv( levels.front()->grid );
    }

    MiniMonT::MiniMonRecord( 1, "smoothfinal" );
}


int main( int argc, char* argv[] ) {

    MiniMonT::MiniMonInit();
    MiniMonT::MiniMonRecord( 0, "main" );

    MiniMonT::MiniMonRecord( 0, "dash::init" );
    dash::init(&argc, &argv);
    auto id= dash::myid();
    MiniMonT::MiniMonRecord( 1, "dash::init" );

    MiniMonT::MiniMonRecord( 0, "setup" );

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

    uint32_t howmanylevels= 5;

    if ( argc > 1 ) {

        if ( 0 == strncmp( "-h", argv[1], 2  ) ||
                0 == strncmp( "--help", argv[1], 6 ) ) {

            if ( 0 == dash::myid() ) {

                cout << "call me as [mpirun] '" << argv[0] << "' [-h|--help] [number-of-levels=5]" << endl;
            }
            exit(0);

        } else {

            /* otherwise interpret as number of grid levels to employ */
            howmanylevels= atoi( argv[1] );
        }
    }

    assert( howmanylevels > 2 );
    assert( howmanylevels <= 16 ); /* please adapt if you really want to go so high */


    vector<Level*> levels;
    levels.reserve( howmanylevels );

    resolutionForCSVd= ( 1<<5 ) * factor_z;
    resolutionForCSVh= ( 1<<5 ) * factor_y;
    resolutionForCSVw= ( 1<<5 ) * factor_x;

    if ( 0 == dash::myid() ) {

        cout << "run '" << argv[0] << "' with " << dash::Team::All().size() << " units "
            "for with grids from " <<
            factor_z << "x" <<
            factor_y << "x" <<
            factor_x <<
            " to " <<
            (1<<(howmanylevels))*factor_z << "x" <<
            (1<<(howmanylevels))*factor_y << "x" <<
            (1<<(howmanylevels))*factor_x <<
            endl << endl;
    }

    /* create all grid levels, starting with the finest and ending with 2x2,
    The finest level is outside the loop because it is always done by dash::Team::All() */

    if ( 0 == dash::myid() ) {
        cout << "finest level is " <<
            (1<<(howmanylevels))*factor_z << "x" <<
            (1<<(howmanylevels))*factor_y << "x" <<
            (1<<(howmanylevels))*factor_x <<
            " distributed over " << dash::Team::All().size() << " units " << endl;
    }

    levels.push_back( new Level(
        (1<<(howmanylevels))*factor_z ,
        (1<<(howmanylevels))*factor_y ,
        (1<<(howmanylevels))*factor_x ,
        dash::Team::All(), teamspec ) );

    initboundary( *levels.back() );
    sanitycheck( levels.back()->grid );

    dash::barrier();

    for ( uint32_t l= 1; l < howmanylevels; l++ ) {

        if ( 0 == dash::myid() ) {
            cout << "compute level " << l << " is " <<
                (1<<(howmanylevels-l))*factor_z << "x" <<
                (1<<(howmanylevels-l))*factor_y << "x" <<
                (1<<(howmanylevels-l))*factor_x <<
                " distributed over " << dash::Team::All().size() << " units " << endl;
        }

        /* do not try to allocate >= 8GB per core -- try to prevent myself
        from running too big a simulation on my laptop */
        assert( (1<<(howmanylevels-l))*factor_z *
            (1<<(howmanylevels-l))*factor_y *
            (1<<(howmanylevels-l))*factor_x < dash::Team::All().size() * (1<<27) );

        levels.push_back(
            new Level( (1<<(howmanylevels-l))*factor_z ,
                        (1<<(howmanylevels-l))*factor_y ,
                        (1<<(howmanylevels-l))*factor_x ,
            dash::Team::All(), teamspec ) );

        initboundary( *levels.back() );

        dash::barrier();
    }

    /* here all units and all teams meet again, those that were active for the coarsest
    levels and those that were dormant */
    dash::Team::All().barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here
    but we do it for demonstration in the graphical output */
    initgrid( levels.front()->grid );

    markunits( levels[0]->grid );

    //writeToCsv( levels.front()->grid );

    dash::Team::All().barrier();

    MiniMonT::MiniMonRecord( 1, "setup" );

    v_cycle( levels.begin(), levels.end(), 10, 0.0001 );
    dash::Team::All().barrier();

    smoothen_final( levels, 0.001 );
    dash::Team::All().barrier();

    writeToCsv( levels.front()->grid );

    MiniMonT::MiniMonRecord( 0, "dash::finalize" );
    dash::finalize();
    MiniMonT::MiniMonRecord( 1, "dash::finalize" );

    MiniMonT::MiniMonRecord( 1, "main" );
    MiniMonT::MiniMonPrint( id, dash::Team::All().size() );

    return 0;
}
