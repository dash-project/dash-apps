#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cstdio>
#include <utility>

#include <libdash.h>
//#include <dash/HaloMatrixWrapper.h>

#include "allreduce.h"
#include "minimonitoring.h"

#define WITHCSVOUTPUT 1
#ifdef WITHCSVOUTPUT

uint32_t filenumber= 0;

#endif /* WITHCSVOUTPUT */


/* TODOs

- add clean version of the code:
    - without asserts
    - without MiniMon
    - with simple loops, no optimization for contiguous lines, etc.

*/


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
using StencilT = dash::Stencil<3>;
using StencilSpecT = dash::StencilSpec<3,6>;
using CycleSpecT = dash::CycleSpec<3>;
using HaloMatrixWrapperT = dash::HaloMatrixWrapper<MatrixT,StencilSpecT>;


const StencilSpecT stencil_spec({
    StencilT(-1, 0, 0), StencilT( 1, 0, 0),
    StencilT( 0,-1, 0), StencilT( 0, 1, 0),
    StencilT( 0, 0,-1), StencilT( 0, 0, 1)});

const CycleSpecT cycle_spec(
    dash::Cycle::FIXED,
    dash::Cycle::FIXED,
    dash::Cycle::FIXED );

struct Level {

    /* now with double-buffering. oldgrid and oldhalo should only be read,
    newgrid should only be written. newhalo is only there to keep newgrid's halo
    before both are swapped in swap() */

    MatrixT* oldgrid;
    HaloMatrixWrapperT* oldhalo;

    MatrixT* newgrid;
    HaloMatrixWrapperT* newhalo;

    Level( size_t d, size_t h, size_t w, dash::Team& team, TeamSpecT teamspec ) {

        oldgrid= new MatrixT( dash::SizeSpec<3>( d, h, w ),
            dash::DistributionSpec<3>( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ),
            team, teamspec );
        newgrid= new MatrixT( dash::SizeSpec<3>( d, h, w ),
            dash::DistributionSpec<3>( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ),
            team, teamspec );
        oldhalo= new HaloMatrixWrapperT( *oldgrid, stencil_spec, cycle_spec );
        newhalo= new HaloMatrixWrapperT( *newgrid, stencil_spec, cycle_spec );
    }

    Level() = delete;

    ~Level() {
/*
        delete oldgrid;
        delete newgrid;
        delete oldhalo;
        delete newhalo;
*/
    }

    /** swap grid and halos for the double buffering scheme */
    void swap() {
// smooth_inner
        std::swap<MatrixT*>( oldgrid, newgrid );
        std::swap<HaloMatrixWrapperT*>( oldhalo, newhalo );
    }


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
void writeToCsv( const MatrixT* grid ) {

/* TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

Here is still a slight shift in the output when the actual grid is larger than the output grid!

*/

#ifdef WITHCSVOUTPUT

    grid->barrier();

    std::array< long int, 3 > corner= grid->pattern().global( {0,0,0} );

    size_t d= grid->extent(0);
    size_t h= grid->extent(1);
    size_t w= grid->extent(2);

    size_t dl= grid->local.extent(0);
    size_t hl= grid->local.extent(1);
    size_t wl= grid->local.extent(2);

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
                    (double) grid->local[ muld*z/divd ][ mulh*y/divh ][ mulw*x/divw ] << "\n";
            }
        }
    }

    csvfile.close();
    grid->barrier();
#endif /* WITHCSVOUTPUT */
}


/* write out the grid in its actual size */
void writeToCsvFullGrid( const MatrixT* grid ) {

#ifdef WITHCSVOUTPUT

    grid->barrier();

    std::array< long int, 3 > corner= grid->pattern().global( {0,0,0} );

    size_t d= grid->extent(0);
    size_t h= grid->extent(1);
    size_t w= grid->extent(2);

    size_t dl= grid->local.extent(0);
    size_t hl= grid->local.extent(1);
    size_t wl= grid->local.extent(2);

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
                    (double) grid->local[ z ][ y ][ x ] << "\n";
            }
        }
    }

    csvfile.close();
    grid->barrier();
#endif /* WITHCSVOUTPUT */
}


void sanitycheck( const MatrixT* grid  ) {

    /* check if the sum of the local extents of the matrix blocks sum up to
    the global extents, abort otherwise */
    dash::Array<size_t> sums( dash::Team::All().size(), dash::BLOCKED );
    sums.local[0]= grid->local.extent(0) * grid->local.extent(1) * grid->local.extent(2);
    sums.barrier();

    if ( 0 != dash::myid() ) {

        dash::transform<size_t>(
            sums.lbegin(), sums.lend(), // first source
            sums.begin(), // second source
            sums.begin(), // destination
            dash::plus<size_t>() );
    }
    sums.barrier();

    if ( ( (size_t) sums[0] ) != grid->extent(0) * grid->extent(1) * grid->extent(2) ) {

        if ( 0 == dash::myid() ) {
            cout << "ERROR: size mismatch: global size is " <<
                grid->extent(0) * grid->extent(1) * grid->extent(2) << " == " <<
                grid->extent(0) << "x" << grid->extent(1) << "x" << grid->extent(2) <<
                " but sum of local sizes is " << (size_t) sums[0] << std::flush << endl;
        }

        dash::finalize();
        exit(0);
    }
}


void initgrid( MatrixT* grid ) {

    /* not strictly necessary but it also avoids NAN values */
    dash::fill( grid->begin(), grid->end(), 5.0 );

    grid->barrier();
}


/* simple boundary settings where 3 sides are hot and 3 are cold */
void initboundary_3hot3cold( Level& level ) {

    double gw= level.oldgrid->extent(2);
    double gh= level.oldgrid->extent(1);
    double gd= level.oldgrid->extent(0);

    /* with this way of setting the boundary conditions on every level separatel (in
    contrast to scaling it down step by step) one needs to make sure that the boundary
    values are continuous (kind of). Otherwise the coarsened version of the boundary
    values will have jumps in different places (instead of jumps being smoothed out
    by coarsening the boundary).
    Well, we should really think about another way to init the boundary values! */
    auto lambda= [gd,gh,gw]( const std::array<dash::default_index_t,3>& coords ) {

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
    };

    level.oldhalo->set_fixed_halos( lambda );
    level.newhalo->set_fixed_halos( lambda );

}


/* alternative boundary settings, the top and bottom planes have a hot circle
in the middle */
// void initboundary_2circles( Level& level ) {
void initboundary( Level& level ) {

    double gw= level.oldgrid->extent(2);
    double gh= level.oldgrid->extent(1);
    double gd= level.oldgrid->extent(0);

    /* This way of setting boundaries uses subsampling on the top and bottom
    planes to determine the border values. This is another logical way that
    may be convenient sometimes. It guarantees that the boundary values on all
    the levels match.
    All other sides are constant at 0.0 degrees. The top an bottom circles are
    hot with 10.0 degrees. */

    auto lambda= [gd,gh,gw]( const std::array<dash::default_index_t,3>& coords ) {

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
    };

    level.oldhalo->set_fixed_halos( lambda );
    level.newhalo->set_fixed_halos( lambda );
}


void markunits( MatrixT* grid ) {

    /* Mark unit bordery by setting the first local rows and columns */

    size_t w= grid->local.extent(2);
    size_t h= grid->local.extent(1);
    size_t d= grid->local.extent(0);

    for ( size_t i = 0; i < d; ++i ) {
        for ( size_t j = 0; j < h; ++j ) {
            grid->local[i][j][0] = 8.0;
        }
    }

    for ( size_t i = 0; i < d; ++i ) {
        for ( size_t k = 0; k < w; ++k ) {
            grid->local[i][0][k] = 8.0;
        }
    }

    for ( size_t j = 0; j < h; ++j ) {
        for ( size_t k = 0; k < w; ++k ) {
            grid->local[0][j][k] = 8.0;
        }
    }
}


void scaledownboundary( Level& fine, Level& coarse ) {

    assert( coarse.oldgrid->extent(2)*2 == fine.oldgrid->extent(2) );
    assert( coarse.oldgrid->extent(1)*2 == fine.oldgrid->extent(1) );
    assert( coarse.oldgrid->extent(0)*2 == fine.oldgrid->extent(0) );

    size_t wc= coarse.oldgrid->local.extent(2);
    size_t hc= coarse.oldgrid->local.extent(1);
    size_t dc= coarse.oldgrid->local.extent(0);
    size_t wf= fine.oldgrid->local.extent(2);
    size_t hf= fine.oldgrid->local.extent(1);
    size_t df= fine.oldgrid->local.extent(0);

    size_t dmax= coarse.oldgrid->extent(0);
    size_t hmax= coarse.oldgrid->extent(1);
    //size_t wmax= coarse.oldgrid->extent(2);

    auto finehalo= fine.oldhalo;

    auto lambda= [&finehalo,&dmax,&hmax]( const std::array<dash::default_index_t,3>& coord ) {

        auto coordf= coord;
        for( auto& c : coordf ) {
            if ( c > 0 ) c *= 2;
        }

        if ( -1 == coord[0] || dmax == coord[0] ) {

            /* z plane */
            return 0.25 * (
                *finehalo->halo_element_at( { coordf[0], coordf[1]+0, coordf[2]+0 } ) +
                *finehalo->halo_element_at( { coordf[0], coordf[1]+0, coordf[2]+1 } ) +
                *finehalo->halo_element_at( { coordf[0], coordf[1]+1, coordf[2]+0 } ) +
                *finehalo->halo_element_at( { coordf[0], coordf[1]+1, coordf[2]+1 } ) );

        } else if ( -1 == coord[1] || hmax == coord[1] ) {

            /* y plane */
            return 0.25 * (
                *finehalo->halo_element_at( { coordf[0]+0, coordf[1], coordf[2]+0 } ) +
                *finehalo->halo_element_at( { coordf[0]+0, coordf[1], coordf[2]+1 } ) +
                *finehalo->halo_element_at( { coordf[0]+1, coordf[1], coordf[2]+0 } ) +
                *finehalo->halo_element_at( { coordf[0]+1, coordf[1], coordf[2]+1 } ) );

        } else /* if ( -1 == coord[2] || wmax == coord[2] ) */ {

            /* x plane */
            return 0.25 * (
                *finehalo->halo_element_at( { coordf[0]+0, coordf[1]+0, coordf[2] } ) +
                *finehalo->halo_element_at( { coordf[0]+0, coordf[1]+1, coordf[2] } ) +
                *finehalo->halo_element_at( { coordf[0]+1, coordf[1]+0, coordf[2] } ) +
                *finehalo->halo_element_at( { coordf[0]+1, coordf[1]+1, coordf[2] } ) );

        }
    };

    coarse.oldhalo->set_fixed_halos( lambda );
    coarse.newhalo->set_fixed_halos( lambda );
}


void scaledown( MatrixT* finegrid, MatrixT* coarsegrid ) {

    uint64_t param= finegrid->local.extent(0)*finegrid->local.extent(1)*finegrid->local.extent(2);
    MiniMonT::MiniMonRecord( 0, "scaledown", param );

    assert( coarsegrid->extent(2) * 2 == finegrid->extent(2) );
    assert( coarsegrid->extent(1) * 2 == finegrid->extent(1) );
    assert( coarsegrid->extent(0) * 2 == finegrid->extent(0) );

    const std::array< long unsigned int, 3 > extentc= coarsegrid->pattern().local_extents();
    const std::array< long signed int, 3 >   cornerc= coarsegrid->pattern().global( {0,0,0} );
    const std::array< long unsigned int, 3 > extentf= finegrid->pattern().local_extents();
    const std::array< long signed int, 3 >   cornerf= finegrid->pattern().global( {0,0,0} );

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

            constexpr size_t x= 0;
            double* __restrict p_coarse= &coarsegrid->local[z][y][x];

            const double* __restrict p_000= &finegrid->local[2*z+0][2*y+0][2*x];
            const double* __restrict p_010= &finegrid->local[2*z+0][2*y+1][2*x];
            const double* __restrict p_100= &finegrid->local[2*z+1][2*y+0][2*x];
            const double* __restrict p_110= &finegrid->local[2*z+1][2*y+1][2*x];

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


void scaleup( MatrixT* coarsegrid, MatrixT* finegrid ) {

    uint64_t param= coarsegrid->local.extent(0)*coarsegrid->local.extent(1)*coarsegrid->local.extent(2);
    MiniMonT::MiniMonRecord( 0, "scaleup", param );

    assert( coarsegrid->extent(2) * 2 == finegrid->extent(2) );
    assert( coarsegrid->extent(1) * 2 == finegrid->extent(1) );
    assert( coarsegrid->extent(0) * 2 == finegrid->extent(0) );

    const std::array< long unsigned int, 3 > extentc= coarsegrid->pattern().local_extents();
    const std::array< long signed int, 3 >   cornerc= coarsegrid->pattern().global( {0,0,0} );
    const std::array< long unsigned int, 3 > extentf= finegrid->pattern().local_extents();
    const std::array< long signed int, 3 >   cornerf= finegrid->pattern().global( {0,0,0} );

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

            constexpr size_t x= 0;
            const double* __restrict p_coarse= &coarsegrid->local[z][y][x];

            double* __restrict p_000= &finegrid->local[2*z+0][2*y+0][2*x];
            double* __restrict p_010= &finegrid->local[2*z+0][2*y+1][2*x];
            double* __restrict p_100= &finegrid->local[2*z+1][2*y+0][2*x];
            double* __restrict p_110= &finegrid->local[2*z+1][2*y+1][2*x];

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


void transfertofewer( Level& source /* with larger team*/, Level& dest /* with smaller team */ ) {

#if 0
    /* should only be called by the smaller team */
    assert( 0 == dest.oldgrid->team().position() );

cout << "unit " << dash::myid() << " transfertofewer" << endl;

    /* we need to find the coordinates that the local unit needs to receive
    from several other units that are not in this team */

    /* we can safely assume that the source blocks are copied entirely */

    std::array< long int, 3 > corner= dest.oldgrid->pattern().global( {0,0,0} );
    std::array< long unsigned int, 3 > sizes= dest.oldgrid->pattern().local_extents();

cout << "    start coord: " <<
    corner[0] << ", "  << corner[1] << ", " << corner[2] << endl;
cout << "    extents: " <<
        sizes[0] << ", "  << sizes[1] << ", " << sizes[2] << endl;
cout << "    dest local  dist " << dest.oldgrid->lend() - dest.oldgrid->lbegin() << endl;
cout << "    dest global dist " << dest.oldgrid->end() - dest.oldgrid->begin() << endl;
cout << "    src  local  dist " << source.oldgrid->lend() - source.oldgrid->lbegin() << endl;
cout << "    src  global dist " << source.oldgrid->end() - source.oldgrid->begin() << endl;

    /* Can I do this any cleverer than loops over the n-1 non-contiguous
    dimensions and then a dash::copy for the 1 contiguous dimension? */

    double buf[512];

    for ( uint32_t z= 0; z < sizes[0]; z++ ) {
        for ( uint32_t y= 0; y < sizes[1]; y++ ) {

            cout << "copy " << corner[0]+z << "," << corner[1]+y << "," <<corner[2] << " -- " <<
                corner[0]+z << "," << corner[1]+y << "," << corner[2] + sizes[2] << " == " <<
                ((corner[0]+z)*sizes[1]+y)*sizes[2] << " - " << ((corner[0]+z)*sizes[1]+y)*sizes[2]+sizes[2] << endl;

            auto start= source.oldgrid->begin() + ((corner[0]+z)*sizes[1]+y)*sizes[2];

            dash::copy( start, start + sizes[2], &dest.oldgrid->local[z][y][0] );
            //dash::copy( start, start + sizes[2], buf );
            //dash::copy( source.grid.begin()+40, source.grid.begin()+48, buf );
        }
    }
#endif /* 0 */
}


void transfertomore( Level& source /* with smaller team*/, Level& dest /* with larger team */ ) {

}


/* Smoothen the given level from oldgrid+oldhalo to newgrid. Call Level::swap() at the end.

This specialization does not compute the residual to have a much simple code. Should be kept in sync with the following version of smoothen() */
void smoothen( Level& level ) {

    size_t ld= level.oldgrid->local.extent(0);
    size_t lh= level.oldgrid->local.extent(1);
    size_t lw= level.oldgrid->local.extent(2);
    uint64_t param= ld*lh*lw;
    MiniMonT::MiniMonRecord( 0, "smoothen", param );

    double res= 0.0;

    /// relaxation coeff.
    const double c= 1.0;

    // async halo update
    level.oldhalo->update_async();

    MiniMonT::MiniMonRecord( 0, "smooth_inner", param );

    // update inner

    /* the start value for both, the y loop and the x loop is 1 because either there is
    a border area next to the halo -- then the first column or row is covered below in
    the border update -- or there is an outside border -- then the first column or row
    contains the boundary values. */
    for ( size_t z= 1; z < ld-1; z++ ) {
        for ( size_t y= 1; y < lh-1; y++ ) {

            /* this should eventually be done with Alpaka or Kokkos to look
            much nicer but still be fast */

            constexpr size_t x= 1;
            const double* __restrict p_here=  &level.oldgrid->local[z  ][y  ][x  ];
            const double* __restrict p_east=  &level.oldgrid->local[z  ][y  ][x+1];
            const double* __restrict p_west=  &level.oldgrid->local[z  ][y  ][x-1];
            const double* __restrict p_north= &level.oldgrid->local[z  ][y+1][x  ];
            const double* __restrict p_south= &level.oldgrid->local[z  ][y-1][x  ];
            const double* __restrict p_up=    &level.oldgrid->local[z+1][y  ][x  ];
            const double* __restrict p_down=  &level.oldgrid->local[z-1][y  ][x  ];

            double* __restrict p_new=   &level.newgrid->local[z  ][y  ][x  ];

            for ( size_t x= 1; x < lw-1; x++ ) {

                double dtheta= ( *p_east + *p_west + *p_north + *p_south + *p_up + *p_down ) / 6.0 - *p_here ;
                *p_new= *p_here + c * dtheta;
                res= std::max( res, std::fabs( dtheta ) );
                p_here++;
                p_east++;
                p_west++;
                p_north++;
                p_south++;
                p_up++;
                p_down++;
                p_new++;
            }
        }
    }

    MiniMonT::MiniMonRecord( 1, "smooth_inner", param );

    MiniMonT::MiniMonRecord( 0, "smooth_wait", param );

    // wait for async halo update
    level.oldhalo->update_async();

    MiniMonT::MiniMonRecord( 1, "smooth_wait", param );

    /* this barrier is there so that iterations are synchronized across all
    units. Otherwise some overtak others esp. on the very small grids. */
    level.oldgrid->barrier();

    MiniMonT::MiniMonRecord( 0, "smooth_outer", param );

    /// begin pointer of local block, needed because halo border iterator is read-only
    auto gridlocalbegin= level.newgrid->lbegin();

    auto bend = level.oldhalo->bend();
    // update border area
    for( auto it = level.oldhalo->bbegin(); it != bend; ++it ) {

        double dtheta= ( it.value_at(0) + it.value_at(1) +
            it.value_at(2) + it.value_at(3) + it.value_at(4) + it.value_at(5) ) / 6.0 - *it;
        gridlocalbegin[ it.lpos() ]= *it + c * dtheta;
        res= std::max( res, std::fabs( dtheta ) );
    }

    MiniMonT::MiniMonRecord( 1, "smooth_outer", param );


    level.swap();

    MiniMonT::MiniMonRecord( 1, "smoothen", param );
}



/**
Smoothen the given level from oldgrid+oldhalo to newgrid. Call Level::swap() at the end.

The parallel global residual is returned as a return parameter, but only
if it is not NULL because then the expensive parallel reduction is just avoided.
*/
double smoothen( Level& level, Allreduce& res ) {

    size_t ld= level.oldgrid->local.extent(0);
    size_t lh= level.oldgrid->local.extent(1);
    size_t lw= level.oldgrid->local.extent(2);
    uint64_t param= ld*lh*lw;
    MiniMonT::MiniMonRecord( 0, "smoothen", param );

    double localres= 0.0;

    /// relaxation coeff.
    const double c= 1.0;

    // async halo update
    level.oldhalo->update_async();

    MiniMonT::MiniMonRecord( 0, "smooth_inner", param );

    // update inner

    /* the start value for both, the y loop and the x loop is 1 because either there is
    a border area next to the halo -- then the first column or row is covered below in
    the border update -- or there is an outside border -- then the first column or row
    contains the boundary values. */
    for ( size_t z= 1; z < ld-1; z++ ) {
        for ( size_t y= 1; y < lh-1; y++ ) {

            /* this should eventually be done with Alpaka or Kokkos to look
            much nicer but still be fast */

            constexpr size_t x= 1;
            const double* __restrict p_here=  &level.oldgrid->local[z  ][y  ][x  ];
            const double* __restrict p_east=  &level.oldgrid->local[z  ][y  ][x+1];
            const double* __restrict p_west=  &level.oldgrid->local[z  ][y  ][x-1];
            const double* __restrict p_north= &level.oldgrid->local[z  ][y+1][x  ];
            const double* __restrict p_south= &level.oldgrid->local[z  ][y-1][x  ];
            const double* __restrict p_up=    &level.oldgrid->local[z+1][y  ][x  ];
            const double* __restrict p_down=  &level.oldgrid->local[z-1][y  ][x  ];

            double* __restrict p_new=   &level.newgrid->local[z  ][y  ][x  ];

            for ( size_t x= 1; x < lw-1; x++ ) {

                double dtheta= ( *p_east + *p_west + *p_north + *p_south + *p_up + *p_down ) / 6.0 - *p_here ;
                *p_new= *p_here + c * dtheta;
                localres= std::max( localres, std::fabs( dtheta ) );
                p_here++;
                p_east++;
                p_west++;
                p_north++;
                p_south++;
                p_up++;
                p_down++;
                p_new++;
            }
        }
    }

    MiniMonT::MiniMonRecord( 1, "smooth_inner", param );

    MiniMonT::MiniMonRecord( 0, "smooth_wait", param );
    // wait for async halo update
    level.oldhalo->wait();
    MiniMonT::MiniMonRecord( 1, "smooth_wait", param );


    MiniMonT::MiniMonRecord( 0, "smooth_col_bc", param );
    /* unit 0 (of any active team) waits until all local residuals from all
    other active units are in */
    res.collect( level.oldgrid->team() );
    res.asyncbroadcast( level.oldgrid->team() );
    MiniMonT::MiniMonRecord( 1, "smooth_col_bc", param );


    /* the former contains a barrier that keeps the iterations in sync */
    // level.oldgrid->barrier();


    MiniMonT::MiniMonRecord( 0, "smooth_outer", param );

    /// begin pointer of local block, needed because halo border iterator is read-only
    auto gridlocalbegin= level.newgrid->lbegin();

    auto bend = level.oldhalo->bend();
    // update border area
    for( auto it = level.oldhalo->bbegin(); it != bend; ++it ) {

        double dtheta= ( it.value_at(0) + it.value_at(1) +
            it.value_at(2) + it.value_at(3) + it.value_at(4) + it.value_at(5) ) / 6.0 - *it;
        gridlocalbegin[ it.lpos() ]= *it + c * dtheta;
        localres= std::max( localres, std::fabs( dtheta ) );
    }

    MiniMonT::MiniMonRecord( 1, "smooth_outer", param );

    MiniMonT::MiniMonRecord( 0, "smooth_wait_set", param );

    res.waitbroadcast( level.oldgrid->team() );

    /* global residual from former iteration */
    double oldres= res.get();

    res.asyncset( localres, level.oldgrid->team() );

    MiniMonT::MiniMonRecord( 1, "smooth_wait_set", param );

    level.swap();

    MiniMonT::MiniMonRecord( 1, "smoothen", param );

    return oldres;
}


/* recursive version */
void v_cycle( vector<Level*>::iterator it, vector<Level*>::iterator itend,
        uint32_t numiter, double epsilon, Allreduce& res ) {

    if ( 0 == dash::myid() ) {
        cout << "v-cycle on  " <<
                    (*it)->oldgrid->extent(2) << "x" <<
                    (*it)->oldgrid->extent(1) << "x" <<
                    (*it)->oldgrid->extent(0) << endl;
    }

    vector<Level*>::iterator itnext( it );
    itnext++;

    /* reached end of recursion? */
    if ( itend == itnext ) {

        /* smoothen completely  */

        uint32_t j= 0;
        while ( res.get() > epsilon ) {

            /* need global residual for iteration count */
            smoothen( **it, res );

            if ( 0 == dash::myid() ) {
                cout << j << ": smoothen coarsest with residual " << res.get() << endl;
            }
            j++;
        }
        writeToCsv( (*it)->oldgrid );
        return;
    }

    /* stepped on the dummy level? ... which is there to signal that it is not
    the end of the parallel recursion on the coarsest level but a subteam is
    going on to solve the coarser levels but this unit is not in that subteam.
    ... sounds complicated, is complicated, change only if you know what you
    are doing. */
    if ( NULL == *itnext ) {

        /* barrier 'Alice', belongs together with the next barrier 'Bob' below */
        (*it)->oldgrid->team().barrier();

        cout << "all meet again here: I'm passive unit " << dash::myid() << endl;
        return;
    }

    /* stepped on a transfer level? */
    if ( (*it)->oldgrid->team().size() != (*itnext)->oldgrid->team().size() ) {

        /* only the members of the reduced team need to work, all others do siesta. */
        //if ( 0 == (*itnext)->grid.team().position() )
        assert( 0 == (*itnext)->oldgrid->team().position() );
        {

            cout << "transfer to " <<
                (*it)->oldgrid->extent(2) << "x" <<
                (*it)->oldgrid->extent(1) << "x" <<
                (*it)->oldgrid->extent(0) << " with " << (*it)->oldgrid->team().size() << " units "
                " --> " <<
                (*itnext)->oldgrid->extent(2) << "x" <<
                (*itnext)->oldgrid->extent(1) << "x" <<
                (*itnext)->oldgrid->extent(0) << " with " << (*itnext)->oldgrid->team().size() << " units " << endl;

            transfertofewer( **it, **itnext );

            v_cycle( itnext, itend, numiter, epsilon, res );

            cout << "transfer back " <<
            (*itnext)->oldgrid->extent(2) << "x" <<
            (*itnext)->oldgrid->extent(1) << "x" <<
            (*itnext)->oldgrid->extent(0) << " with " << (*itnext)->oldgrid->team().size() << " units "
            " --> " <<
            (*it)->oldgrid->extent(2) << "x" <<
            (*it)->oldgrid->extent(1) << "x" <<
            (*it)->oldgrid->extent(0) << " with " << (*it)->oldgrid->team().size() << " units " <<  endl;

            transfertomore( **itnext, **it );
        }

        /* barrier 'Bob', belongs together with the previous barrier 'Alice' above */
        (*it)->oldgrid->team().barrier();


        cout << "all meet again here: I'm active unit " << dash::myid() << endl;
        return;
    }

    /* normal recursion */

    /* smoothen somewhat with fixed number of iterations */
    for ( uint32_t j= 0; j < numiter; j++ ) {
        smoothen( **it );
    }

    writeToCsv( (*it)->oldgrid );

    /* scale down */
    if ( 0 == dash::myid() ) {
        cout << "scale down " <<
            (*it)->oldgrid->extent(2) << "x" <<
            (*it)->oldgrid->extent(1) << "x" <<
            (*it)->oldgrid->extent(0) <<
            " --> " <<
            (*itnext)->oldgrid->extent(2) << "x" <<
            (*itnext)->oldgrid->extent(1) << "x" <<
            (*itnext)->oldgrid->extent(0) << endl;
    }
    scaledown( (*it)->oldgrid, (*itnext)->oldgrid );
    writeToCsv( (*itnext)->oldgrid );

    /* recurse  */
    v_cycle( itnext, itend, numiter, epsilon, res );

    /* scale up */
    if ( 0 == dash::myid() ) {
        cout << "scale up " <<
            (*itnext)->oldgrid->extent(2) << "x" <<
            (*itnext)->oldgrid->extent(1) << "x" <<
            (*itnext)->oldgrid->extent(0) <<
            " --> " <<
            (*it)->oldgrid->extent(2) << "x" <<
            (*it)->oldgrid->extent(1) << "x" <<
            (*it)->oldgrid->extent(0) << endl;
    }
    scaleup( (*itnext)->oldgrid, (*it)->oldgrid );
    writeToCsv( (*it)->oldgrid );

    /* smoothen somewhat with fixed number of iterations */
    for ( uint32_t j= 0; j < numiter; j++ ) {
        smoothen( **it );
    }
    writeToCsv( (*it)->oldgrid );
}


void smoothen_final( Level* level, double epsilon, Allreduce& res ) {

    MiniMonT::MiniMonRecord( 0, "smoothfinal" );

    uint32_t j= 0;
    while ( res.get() > epsilon ) {

        smoothen( *level, res );

        if ( 0 == dash::myid() ) {
            cout << j << ": smoothen finest with residual " << res.get() << endl;
        }
        j++;
    }

    MiniMonT::MiniMonRecord( 1, "smoothfinal" );
}


void do_multigrid_iteration( uint32_t howmanylevels ) {

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

    vector<Level*> levels;
    levels.reserve( howmanylevels );

    resolutionForCSVd= ( 1<<6 ) * factor_z/factor_max;
    resolutionForCSVh= ( 1<<6 ) * factor_y/factor_max;
    resolutionForCSVw= ( 1<<6 ) * factor_x/factor_max;

    if ( 0 == dash::myid() ) {

        cout << "run multigrid iteration with " << dash::Team::All().size() << " units "
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

    /* only do initgrid on the finest level, use caledownboundary for all others */
    initboundary( *levels.back() );
    sanitycheck( levels.back()->oldgrid );

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

        Level& previouslevel= *levels.back();

        levels.push_back(
            new Level( (1<<(howmanylevels-l))*factor_z ,
                        (1<<(howmanylevels-l))*factor_y ,
                        (1<<(howmanylevels-l))*factor_x ,
            dash::Team::All(), teamspec ) );

        /* scaledown boundary instead of initializing it from the same
        procedure, because this is very prone to subtle mistakes which
        makes the entire multigrid algorithm misbehave. */
        scaledownboundary( previouslevel, *levels.back() );

        dash::barrier();
    }

    /* here all units and all teams meet again, those that were active for the coarsest
    levels and those that were dormant */
    dash::Team::All().barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here
    but we do it for demonstration in the graphical output */
    initgrid( levels.front()->oldgrid );
    markunits( levels.front()->oldgrid );
    writeToCsv( levels.front()->oldgrid );

    dash::Team::All().barrier();

    Allreduce res( dash::Team::All() );
    res.reset( dash::Team::All() );

    MiniMonT::MiniMonRecord( 1, "setup" );

    v_cycle( levels.begin(), levels.end(), 10, 0.001, res );
    dash::Team::All().barrier();
    v_cycle( levels.begin(), levels.end(), 10, 0.0001, res );
    dash::Team::All().barrier();

    res.reset( dash::Team::All() );
    smoothen_final( levels.front(), 0.001, res );
    for ( int i= 0; i < 5; ++i ) {
        writeToCsv( levels.front()->oldgrid );
    }

    dash::Team::All().barrier();
}



void do_multigrid_elastic( uint32_t howmanylevels ) {

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

    vector<Level*> levels;
    levels.reserve( howmanylevels );

    resolutionForCSVd= ( 1<<6 ) * factor_z/factor_max;
    resolutionForCSVh= ( 1<<6 ) * factor_y/factor_max;
    resolutionForCSVw= ( 1<<6 ) * factor_x/factor_max;

    if ( 0 == dash::myid() ) {

        cout << "run elastic multigrid iteration with " << dash::Team::All().size() << " units "
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

    /* only do initgrid on the finest level, use caledownboundary for all others */
    initboundary( *levels.back() );
    sanitycheck( levels.back()->oldgrid );

    dash::barrier();

    for ( uint32_t l= 1; l < howmanylevels; l++ ) {

        dash::Team& previousteam= levels.back()->oldgrid->team();
        dash::Team& currentteam= ( 4 == l ) ? previousteam.split(4) : previousteam;
        TeamSpecT localteamspec( currentteam.size(), 1, 1 );
        localteamspec.balance_extents();

        if ( 0 == currentteam.position() ) {

            if ( previousteam.size() != currentteam.size() ) {

                /* the team working on the following grid layers has just
                been reduced. Therefore, we add an additional grid with the
                same size as the previous one but for the reduced team. Then,
                copying the data from the domain of the larger team to the
                domain of the smaller team is easy. */

                if ( 0 == currentteam.myid() ) {
                    cout << "transfer level " << l-1 << " is " <<
                        (1<<(howmanylevels-l+1))*factor_z << "x" <<
                        (1<<(howmanylevels-l+1))*factor_y << "x" <<
                        (1<<(howmanylevels-l+1))*factor_x <<
                        " distributed over " <<
                        localteamspec.num_units(0) << "x" <<
                        localteamspec.num_units(1) << "x" <<
                        localteamspec.num_units(2) << " units ()" << endl;
                }

                levels.push_back(
                    new Level( (1<<(howmanylevels-l+1))*factor_z ,
                                (1<<(howmanylevels-l+1))*factor_y ,
                                (1<<(howmanylevels-l+1))*factor_x ,
                    currentteam, localteamspec ) );
            }

            /*
            cout << "working unit " << dash::myid() << " / " << currentteam.myid() << " in subteam at position " << currentteam.position() << endl;
            */

            if ( 0 == currentteam.myid() ) {
                cout << "compute level " << l << " is " <<
                    (1<<(howmanylevels-l))*factor_z << "x" <<
                    (1<<(howmanylevels-l))*factor_y << "x" <<
                    (1<<(howmanylevels-l))*factor_x <<
                    " distributed over " <<
                    localteamspec.num_units(0) << "x" <<
                    localteamspec.num_units(1) << "x" <<
                    localteamspec.num_units(2) << " units ()" << endl;
            }

            /* do not try to allocate >= 8GB per core -- try to prevent myself
            from running too big a simulation on my laptop */
            assert( (1<<(howmanylevels-l))*factor_z *
                (1<<(howmanylevels-l))*factor_y *
                (1<<(howmanylevels-l))*factor_x < currentteam.size() * (1<<27) );

            levels.push_back(
                new Level( (1<<(howmanylevels-l))*factor_z ,
                            (1<<(howmanylevels-l))*factor_y ,
                            (1<<(howmanylevels-l))*factor_x ,
                currentteam, localteamspec ) );

            currentteam.barrier();

            /*
            cout << "unit " << dash::myid() << " : " <<
                levels.back()->grid.local.extent(1) << " x " <<
                levels.back()->grid.local.extent(0) << endl;
            dash::barrier();
            */

            /* should use scaledownboundary() but it is more complicated here
            ... let's find out which is better in the end */
            initboundary( *levels.back() );

        } else {

            /*
            cout << "waiting unit " << dash::myid() << " / " << currentteam.myid() << " in subteam at position " << currentteam.position() << endl;
            */

            /* this is a passive unit not takin part in the subteam that
            handles the coarser grids. insert a dummy entry in the vector
            of levels to signal that this is not the coarsest level globally. */
            levels.push_back( NULL );

            break;
        }
    }

    /* here all units and all teams meet again, those that were active for the coarsest
    levels and those that were dormant */
    dash::Team::All().barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here
    but we do it for demonstration in the graphical output */
    initgrid( levels.front()->oldgrid );
    markunits( levels.front()->oldgrid );
    writeToCsv( levels.front()->oldgrid );

    dash::Team::All().barrier();

    Allreduce res( dash::Team::All() );
    res.reset( dash::Team::All() );

    MiniMonT::MiniMonRecord( 1, "setup" );

    v_cycle( levels.begin(), levels.end(), 10, 0.001, res );
    dash::Team::All().barrier();
    v_cycle( levels.begin(), levels.end(), 10, 0.0001, res );
    dash::Team::All().barrier();

    res.reset( dash::Team::All() );
    smoothen_final( levels.front(), 0.001, res );
    for ( int i= 0; i < 5; ++i ) {
        writeToCsv( levels.front()->oldgrid );
    }

    dash::Team::All().barrier();
}


void do_flat_iteration( uint32_t howmanylevels ) {


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

    vector<Level*> levels;
    levels.reserve( howmanylevels );

    resolutionForCSVd= ( 1<<6 ) * factor_z/factor_max;
    resolutionForCSVh= ( 1<<6 ) * factor_y/factor_max;
    resolutionForCSVw= ( 1<<6 ) * factor_x/factor_max;

    if ( 0 == dash::myid() ) {

        cout << "run flat iteration with " << dash::Team::All().size() << " units "
            "for grids of " <<
            (1<<(howmanylevels))*factor_z << "x" <<
            (1<<(howmanylevels))*factor_y << "x" <<
            (1<<(howmanylevels))*factor_x <<
            endl << endl;
    }

    Level* level= new Level(
        (1<<(howmanylevels))*factor_z ,
        (1<<(howmanylevels))*factor_y ,
        (1<<(howmanylevels))*factor_x ,
        dash::Team::All(), teamspec );

    dash::barrier();

    initboundary( *level );
    initgrid( level->oldgrid );
    markunits( level->oldgrid );
    writeToCsv( level->oldgrid );

    dash::barrier();

    //sanitycheck( level->oldgrid );

    //dash::barrier();

    MiniMonT::MiniMonRecord( 0, "smoothflatfixed" );

    uint32_t j= 0;
    while ( j < 200 ) {

        smoothen( *level );

        if ( 0 == dash::myid() ) {
            cout << j << ": smoothen grid without residual " << endl;
        }
        j++;
        //writeToCsv( level->oldgrid );
    }

    MiniMonT::MiniMonRecord( 1, "smoothflatfixed" );

    writeToCsv( level->oldgrid );

    MiniMonT::MiniMonRecord( 0, "smoothflatresidual" );

    Allreduce res( dash::Team::All() );
    res.reset( dash::Team::All() );

    double epsilon= 0.001;
    while ( res.get() > epsilon && j < 400 ) {

        smoothen( *level, res );

        if ( 0 == dash::myid() ) {
            cout << j << ": smoothen grid with residual " << res.get() << endl;
        }
        j++;
        //writeToCsv( level->oldgrid );
    }

    MiniMonT::MiniMonRecord( 1, "smoothflatresidual" );

    writeToCsv( level->oldgrid );

    delete level;
    level= NULL;
}


int main( int argc, char* argv[] ) {

    MiniMonT::MiniMonInit();
    MiniMonT::MiniMonRecord( 0, "main" );

    MiniMonT::MiniMonRecord( 0, "dash::init" );
    dash::init(&argc, &argv);
    auto id= dash::myid();
    MiniMonT::MiniMonRecord( 1, "dash::init" );

    MiniMonT::MiniMonRecord( 0, "setup" );


    bool do_flatsolver= false;
    bool do_elastic= false;
    uint32_t howmanylevels= 5;


    for ( int a= 1; a < argc; a++ ) {

        if ( 0 == strncmp( "-h", argv[a], 2  ) ||
                0 == strncmp( "--help", argv[a], 6 ) ) {

            if ( 0 == dash::myid() ) {

                cout << "call me as [mpirun] '" << argv[0] << "' [-h|--help] [number-of-levels=5]" << endl;
            }
            exit(0);

        } else if ( 0 == strncmp( "-f", argv[a], 2  ) ||
                0 == strncmp( "--flat", argv[a], 7 )) {

            do_flatsolver= true;
            if ( 0 == dash::myid() ) {

                cout << "do flat iteration instead of multigrid" << endl;
            }

        } else if ( 0 == strncmp( "-e", argv[a], 2  ) ||
                0 == strncmp( "--elastic", argv[a], 7 )) {

            do_elastic= true;
            if ( 0 == dash::myid() ) {

                cout << "do multigrid iteration with changing number of units per grid" << endl;
            }

        } else {

            /* otherwise interpret as number of grid levels to employ */
            howmanylevels= atoi( argv[a] );
            if ( 0 == dash::myid() ) {
                cout << "using " << howmanylevels << " levels, " <<
                (1<<howmanylevels) << "^3" << " per unit" << endl;
            }
        }
    }

    assert( howmanylevels > 2 );
    assert( howmanylevels <= 16 ); /* please adapt if you really want to go so high */

    if ( do_flatsolver ) {

        do_flat_iteration( howmanylevels );

    } else if ( do_elastic ) {

        do_multigrid_elastic( howmanylevels );

    } else {

        do_multigrid_iteration( howmanylevels );
    }

    MiniMonT::MiniMonRecord( 0, "dash::finalize" );
    dash::finalize();
    MiniMonT::MiniMonRecord( 1, "dash::finalize" );

    MiniMonT::MiniMonRecord( 1, "main" );
    MiniMonT::MiniMonPrint( id, dash::Team::All().size() );

    return 0;
}
