#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cstdio>
#include <utility>
#include <math.h>

#include "allreduce.h"
#include "minimonitoring.h"

#ifdef WITHCSVOUTPUT

/* make it a shared value to keep it in sync over all units
even in elastic mode when some units take part in some iterations
and some do not. Having everyone ++ it individually won't work anymore. */
dash::Shared<uint32_t>* filenumber;


#endif /* WITHCSVOUTPUT */

/* TODOs

- add clean version of the code:
    - without asserts
    - with simple loops, no optimization for contiguous lines, etc.

*/

MiniMon minimon;

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
using StencilSpecT = dash::StencilSpec<3,14>;
using CycleSpecT = dash::CycleSpec<3>;

/* for the smoothing operation, only the 6-point stencil is needed. 
However, the prolongation operation also needs the */
constexpr StencilSpecT stencil_spec({
    StencilT(-1, 0, 0), StencilT( 1, 0, 0),
    StencilT( 0,-1, 0), StencilT( 0, 1, 0),
    StencilT( 0, 0,-1), StencilT( 0, 0, 1),

    StencilT(-1,-1,-1), StencilT( 1,-1,-1),
    StencilT(-1,-1, 1), StencilT( 1,-1, 1),
    StencilT(-1, 1,-1), StencilT( 1, 1,-1),
    StencilT(-1, 1, 1), StencilT( 1, 1, 1)});

constexpr CycleSpecT cycle_spec(
    dash::Cycle::FIXED,
    dash::Cycle::FIXED,
    dash::Cycle::FIXED );

struct Level {

public:
  using SizeSpecT = dash::SizeSpec<3>;
  using DistSpecT = dash::DistributionSpec<3>;
  using HaloT = dash::HaloMatrixWrapper<MatrixT,StencilSpecT>;


    /* now with double-buffering. src_grid and src_halo should only be read,
    newgrid should only be written. dst_grid and dst_halo are only there to keep the other ones
    before both are swapped in swap() */

public:
    MatrixT* src_grid;
    MatrixT* dst_grid;
    MatrixT* rhs_grid; /* right hand side, doesn't need a halo */
    HaloT* src_halo;
    HaloT* dst_halo;

    /* this are the values of the 7 non-zero matrix values -- only 4 different values, though,
    because the matrix is symmetric */
    double acenter, ax, ay, az;
    /* this is the factor for the right hand side, which is 0 at the finest grid. */
    double ff;
    /* Diagonal element of matrix M, which is the inverse of the diagonal of matrix A. 
    This factor multiplies the defect in $ f - Au $. */
    double m;

    /* sz, sy, sx are the dimensions in meters of the grid excluding the boundary regions */
    double sz, sy, sx;

    /***
    lz, ly, lx are the dimensions in meters of the grid including the boundary regions,
    nz, ny, nx are th number of inner grid points per dimension, excluding the boundary regions,
    therefore, lz,ly,lx are discretized into (nz+2)*(ny+2)*(nx+2) grid points
    */
    Level( double lz, double ly, double lx,
           size_t nz, size_t ny, size_t nx,
           dash::Team& team, TeamSpecT teamspec ) :
            _grid_1( SizeSpecT( nz, ny, nx ), DistSpecT( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
            _grid_2( SizeSpecT( nz, ny, nx ), DistSpecT( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
            _rhs_grid( SizeSpecT( nz, ny, nx ), DistSpecT( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
        _halo_grid_1( _grid_1, stencil_spec, cycle_spec ),
        _halo_grid_2( _grid_2, stencil_spec, cycle_spec ),
        src_grid(&_grid_1), dst_grid(&_grid_2), rhs_grid(&_rhs_grid),
        src_halo(&_halo_grid_1), dst_halo(&_halo_grid_2) {

        assert( 1 < nz );
        assert( 1 < ny );
        assert( 1 < nx );

        sz= lz;
        sy= ly;
        sx= lx;

        double hz= lz/(nz+1);
        double hy= ly/(ny+1);
        double hx= lx/(nx+1);

        /* This is the original setting for the linear system. */

        ax= -1.0/hx/hx;
        ay= -1.0/hy/hy;
        az= -1.0/hz/hz;
        acenter= -2.0*(ax+ay+az);
        m= 1.0 / acenter;

        ff= 1.0; /* factor for right-hand-side */

        for ( uint32_t a= 0; a < team.size(); a++ ) {
            if ( a == dash::myid() ) {

                if ( 0 == a ) {
                    cout << "Level " <<
                        "dim. " << lz << "m*" << ly << "m*" << lz << "m " <<
                        "in grid of " << nz << "*" << ny << "*" << nx <<
                        " h_= " << hz << "," << hy << "," << hx <<
                        " with team of " << team.size() << 
                        " --> a_= " << acenter << "," << ax << "," << ay << "," << az << 
                        " , m= " << m << " , ff= " << ff <<endl;
                }
            }

            team.barrier();
        }
    }

    /***
    Alternative version of the constructor that takes the parent Level as the first argument.
    From this, it can get the original physical dimensions lz, ly, lx and the original
    grid distances hy, hy, hx.
    nz, ny, nx are th number of inner grid points per dimension, excluding the boundary regions
    */
    Level( const Level& parent,
           size_t nz, size_t ny, size_t nx,
           dash::Team& team, TeamSpecT teamspec ) :
            _grid_1( SizeSpecT( nz, ny, nx ), DistSpecT( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
            _grid_2( SizeSpecT( nz, ny, nx ), DistSpecT( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
            _rhs_grid( SizeSpecT( nz, ny, nx ), DistSpecT( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
        _halo_grid_1( _grid_1, stencil_spec, cycle_spec ),
        _halo_grid_2( _grid_2, stencil_spec, cycle_spec ),
        src_grid(&_grid_1), dst_grid(&_grid_2), rhs_grid(&_rhs_grid),
        src_halo(&_halo_grid_1), dst_halo(&_halo_grid_2) {

        assert( 1 < nz );
        assert( 1 < ny );
        assert( 1 < nx );

        sz= parent.sz;
        sy= parent.sy;
        sx= parent.sx;

        ax= parent.ax;
        ay= parent.ay;
        az= parent.az;
        acenter= parent.acenter;
        ff= parent.ff;
        m= parent.m;

        for ( uint32_t a= 0; a < team.size(); a++ ) {
            if ( a == dash::myid() ) {

                if ( 0 == a ) {
                    cout << "Level with a parent level " <<
                        "in grid of " << nz << "*" << ny << "*" << nx <<
                        " with team of " << team.size() << 
                        " --> a_= " << acenter << "," << ax << "," << ay << "," << az <<
                        " , m= " << m << " , ff= " << ff << endl;
                }
            }

            team.barrier();
        }
    }

    Level() = delete;

    /** swap grid and halos for the double buffering scheme */
    void swap() {

        std::swap( src_halo, dst_halo );
        std::swap( src_grid, dst_grid );
    }

    /* to be called by the master unit alone */
    void printout() {

        src_grid->barrier();

        auto n= src_grid->extents();

        for ( size_t z= 0; z < n[0]; ++z ) {
            cout << z << ")" << endl;
            for ( size_t y= 0; y < n[1]; ++y ) {
                for ( size_t x= 0; x < n[2]; ++x ) {

                    cout << std::setprecision(5) << (double) (*src_grid)[z][y][x] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }

    }

    /* to be called by the master unit alone */
    void printout_dst() {

        dst_grid->barrier();

        auto n= dst_grid->extents();

        for ( size_t z= 0; z < n[0]; ++z ) {
            cout << z << ")" << endl;
            for ( size_t y= 0; y < n[1]; ++y ) {
                for ( size_t x= 0; x < n[2]; ++x ) {

                    cout << std::setprecision(5) << (double) (*dst_grid)[z][y][x] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }

    }

    /* to be called by the master unit alone */
    void printout_rhs() {

        rhs_grid->barrier();

        auto n= rhs_grid->extents();

        for ( size_t z= 0; z < n[0]; ++z ) {
            cout << z << ")" << endl;
            for ( size_t y= 0; y < n[1]; ++y ) {
                for ( size_t x= 0; x < n[2]; ++x ) {

                    cout << std::setprecision(5) << (double) (*rhs_grid)[z][y][x] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }

    }

    /* to be called by all units together */
    void printout_halo() {

        std::array< size_t, 3 > n= src_grid->local.extents();
        std::array< long int, 3 > corner= src_grid->pattern().global( {0,0,0} );

        if ( 0 == dash::myid() ) {
            cout << " === z= -1 === " << endl; 
        }

        for ( uint32_t unit= 0; unit < dash::Team::All().size(); ++unit ) {

            src_grid->barrier();
            if ( dash::myid() == unit ) {

                cout << endl << "unit " << unit << ": " << n[1] << "," << n[2] << endl; 
                int z= -1;
                for ( int y= -1; y <= (int) n[1]; ++y ) {
                    for ( int x= -1; x <= (int) n[2]; ++x ) {

                        cout << std::setprecision(3) << 
                            (double) *src_halo->halo_element_at( { corner[0]+z, corner[1]+y, corner[2]+x } ) << " ";

                    }
                    cout << endl;
                }
                cout << endl;
            }
            src_grid->barrier();
        }

        if ( 0 == dash::myid() ) {
            cout << " === z= max === " << endl; 
        }

        for ( uint32_t unit= 0; unit < dash::Team::All().size(); ++unit ) {

            src_grid->barrier();
            if ( dash::myid() == unit ) {

                cout << endl << "unit " << unit << ": " << n[1] << "," << n[2] << endl; 
                int z= n[0];
                for ( int y= -1; y <= (int) n[1]; ++y ) {
                    for ( int x= -1; x <= (int) n[2]; ++x ) {

                        cout << std::setprecision(3) << 
                            (double) *src_halo->halo_element_at( { corner[0]+z, corner[1]+y, corner[2]+x } ) << " ";

                    }
                    cout << endl;
                }
                cout << endl;
            }
            src_grid->barrier();
        }

        if ( 0 == dash::myid() ) {
            cout << " === y= -1 === " << endl; 
        }

        for ( uint32_t unit= 0; unit < dash::Team::All().size(); ++unit ) {

            src_grid->barrier();
            if ( dash::myid() == unit ) {

                cout << endl << "unit " << unit << ": " << n[0] << "," << n[2] << endl; 
                int y= -1;
                for ( int z= -1; z <= (int) n[0]; ++z ) {
                    for ( int x= -1; x <= (int) n[2]; ++x ) {

                        cout << std::setprecision(3) << 
                            (double) *src_halo->halo_element_at( { corner[0]+z, corner[1]+y, corner[2]+x } ) << " ";

                    }
                    cout << endl;
                }
                cout << endl;
            }
            src_grid->barrier();
        }

        if ( 0 == dash::myid() ) {
            cout << " === y= max === " << endl; 
        }

        for ( uint32_t unit= 0; unit < dash::Team::All().size(); ++unit ) {

            src_grid->barrier();
            if ( dash::myid() == unit ) {

                cout << endl << "unit " << unit << ": " << n[0] << "," << n[2] << endl; 
                int y= n[1];
                for ( int z= -1; z <= (int) n[0]; ++z ) {
                    for ( int x= -1; x <= (int) n[2]; ++x ) {

                        cout << std::setprecision(3) << 
                            (double) *src_halo->halo_element_at( { corner[0]+z, corner[1]+y, corner[2]+x } ) << " ";

                    }
                    cout << endl;
                }
                cout << endl;
            }
            src_grid->barrier();
        }

        if ( 0 == dash::myid() ) {
            cout << " === x= -1 === " << endl; 
        }

        for ( uint32_t unit= 0; unit < dash::Team::All().size(); ++unit ) {

            src_grid->barrier();
            if ( dash::myid() == unit ) {

                cout << endl << "unit " << unit << ": " << n[0] << "," << n[1] << endl; 
                int x= -1;
                for ( int z= -1; z <= (int) n[0]; ++z ) {
                    for ( int y= -1; y <= (int) n[1]; ++y ) {

                        cout << std::setprecision(3) << 
                            (double) *src_halo->halo_element_at( { corner[0]+z, corner[1]+y, corner[2]+x } ) << " ";

                    }
                    cout << endl;
                }
                cout << endl;
            }
            src_grid->barrier();
        }

        if ( 0 == dash::myid() ) {
            cout << " === x= max === " << endl; 
        }

        for ( uint32_t unit= 0; unit < dash::Team::All().size(); ++unit ) {

            src_grid->barrier();
            if ( dash::myid() == unit ) {

                cout << endl << "unit " << unit << ": " << n[0] << "," << n[1] << endl; 
                int x= n[2];
                for ( int z= -1; z <= (int) n[0]; ++z ) {
                    for ( int y= -1; y <= (int) n[1]; ++y ) {

                        cout << std::setprecision(3) << 
                            (double) *src_halo->halo_element_at( { corner[0]+z, corner[1]+y, corner[2]+x } ) << " ";

                    }
                    cout << endl;
                }
                cout << endl;
            }
            src_grid->barrier();
        }

    }

private:
    MatrixT _grid_1;
    MatrixT _grid_2;
    HaloT _halo_grid_1;
    HaloT _halo_grid_2;
    MatrixT _rhs_grid;
};


/* global resolution for cvs output, should be fixed such that 
paraview gets input of constant dimensions. Any size > 2 should be good
but odd numbers are suggested. */
std::array< long int, 3 > resolution= {33,33,33};

/* helper function for the following write_to_cvs() function */
inline double arbitrary_element( const MatrixT& grid, Level::HaloT& halo, 
        std::array< long int, 3 >& corner, std::array< size_t, 3 >& localdim, 
        int zz, int yy, int xx ) {

    /* is halo value? That is when any coordinate is -1 or dim[.]. 
    TODO replace by halo convenience layer that provides either 
    halo value or grid value for coordinates [-1:dim[.]] inclusively. */
    if ( -1 == zz-corner[0] || -1 == yy-corner[1] || -1 == xx-corner[2] || 
            localdim[0] == zz-corner[0] || localdim[1] == yy-corner[1] || localdim[2] == xx-corner[2] ) {

        /* get halo value */
        return *halo.halo_element_at( {zz,yy,xx} );
    } else {

        /* get grid value */
        return grid.local[zz-corner[0]][yy-corner[1]][xx-corner[2]];
    }
}


/* write out the current state of the given grid as CSV so that Paraview can
read and visualize it as a structured grid. Use a fixed size for the grid as
defined by the previous global variable. So an animations of the
multigrid procedure is possible.

Do not restrict the grid sizes in any way. Lay an output grid over the source grid, 
then interpolate every output grid point from the 8 neighbor source grid points.

Require halos, include boundary conditions in output.

Accept checks whether any point is part of the grid or part of the halo to make this code simpler.

*/
//void writeToCsv_interpolate( const Level& level ) {
void writeToCsv( const Level& level ) {

#ifdef WITHCSVOUTPUT

    using signed_size_t = typename std::make_signed<size_t>::type;

    const MatrixT& grid= *level.src_grid;
    Level::HaloT& halo= *level.src_halo;

    std::ofstream csvfile;
    std::ostringstream num_string;
    num_string << std::setw(5) << std::setfill('0') << (uint32_t) filenumber->get();
    csvfile.open( "image_unit" + std::to_string(grid.team().myid()) +
        ".csv." + num_string.str() );

    std::array< long int, 3 > corner= grid.pattern().global( {0,0,0} );
    std::array< size_t, 3 > dim= grid.extents();
    std::array< size_t, 3 > localdim= grid.local.extents();

    /*
    cout << "unit " << dash::myid() << " localdim " << 
        localdim[0] << "," << localdim[1] << "," << localdim[2] << endl;
    */

    std::array< long int, 3 > corner0= grid.pattern().global( {0,0,0} );
    std::array< long int, 3 > corner1= 
        grid.pattern().global( {(signed_size_t)localdim[0]-1,(signed_size_t)localdim[1]-1,(signed_size_t)localdim[2]-1} );

    /* start and stop are in output grid coordinates */
    std::array< long int, 3 > start;
    std::array< long int, 3 > stop;
    for ( size_t i= 0; i < 3; ++i ) {

        start[i]= corner[i] * resolution[i] / dim[i];
        stop[i]= ( corner[i] + localdim[i] ) * resolution[i] / dim[i];
    }

    /*
    cout << "unit " << dash::myid() << " range " << 
        start[0] << "," << start[1] << "," << start[2] << " - " <<
        stop[0] << "," << stop[1] << "," << stop[2] << endl;
    */

    grid.barrier();
    if ( 0 == dash::myid() ) {
        csvfile << " z-index,y-index,x-index,z-coord,y-coord,x-coord,heat" << "\n";
        filenumber->set( 1 + (uint32_t) filenumber->get()  );
    }
    grid.barrier();

    /* update halo values, cannot do it async here because there is not much else to do. */
    halo.update();

    /* z,y,z are in output grid coordinates */
    for ( int z= start[0]; z < stop[0]; ++z ) {
        for ( int y= start[1]; y < stop[1]; ++y ) {
            for ( int x= start[2]; x < stop[2]; ++x ) {

                /* for every point of the output grid, determine the 8 nearest
                input grid points. They are represented as pos[]+add[] where 
                add= {1,1,1} except for border cases. */

                /* pos and pos_double are in input grid coordinates */
                std::array< long int, 3 > pos;
                pos[0]= z * (dim[0]+1) / (resolution[0]-1) -1;
                pos[1]= y * (dim[1]+1) / (resolution[1]-1) -1;
                pos[2]= x * (dim[2]+1) / (resolution[2]-1) -1;

                /* factorf(ront) and factor_b(ack) with factor_f[.] + factor_b[.] == 1.0 i.e., 
                the linear interpolation per dimension adds to 1.0 always */
                std::array< double, 3 > factor_b;;
                factor_b[0]= ((double) z) * (dim[0]+1.0) / (resolution[0]-1) - 1.0 - pos[0];
                factor_b[1]= ((double) y) * (dim[1]+1.0) / (resolution[1]-1) - 1.0 - pos[1];
                factor_b[2]= ((double) x) * (dim[2]+1.0) / (resolution[2]-1) - 1.0 - pos[2];

                assert( 0.0 <= factor_b[0] && factor_b[0] < 1.0 );
                assert( 0.0 <= factor_b[1] && factor_b[1] < 1.0 );
                assert( 0.0 <= factor_b[2] && factor_b[2] < 1.0 );

                std::array< double, 3 > factor_f;
                factor_f[0]= 1.0 - factor_b[0];
                factor_f[1]= 1.0 - factor_b[1];
                factor_f[2]= 1.0 - factor_b[2];

                std::array< long int, 3 > add;
                add[0]= ( 0 == z || stop[0]-1 == z ) ? 0 : 1;
                add[1]= ( 0 == y || stop[1]-1 == y ) ? 0 : 1;
                add[2]= ( 0 == x || stop[2]-1 == x ) ? 0 : 1;

                double value= 0.0;

                value += factor_f[0]*factor_f[1]*factor_f[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]       , pos[1]       , pos[2]        );
                value += factor_f[0]*factor_f[1]*factor_b[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]       , pos[1]       , pos[2]+add[2] );
                value += factor_f[0]*factor_b[1]*factor_f[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]       , pos[1]+add[1], pos[2]        );
                value += factor_f[0]*factor_b[1]*factor_b[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]       , pos[1]+add[1], pos[2]+add[2] );
                value += factor_b[0]*factor_f[1]*factor_f[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]+add[0], pos[1]       , pos[2]        );
                value += factor_b[0]*factor_f[1]*factor_b[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]+add[0], pos[1]       , pos[2]+add[2] );
                value += factor_b[0]*factor_b[1]*factor_f[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]+add[0], pos[1]+add[1], pos[2]        );
                value += factor_b[0]*factor_b[1]*factor_b[2] * arbitrary_element( grid, halo, corner, localdim, pos[0]+add[0], pos[1]+add[1], pos[2]+add[2] );

                csvfile << 
                    setfill('0') << setw(4) << z << "," <<
                    setfill('0') << setw(4) << y << "," <<
                    setfill('0') << setw(4) << x << "," <<
                    (double) level.sz * z / (resolution[0]-1) << "," <<
                    (double) level.sy * y / (resolution[1]-1) << "," <<
                    (double) level.sx * x / (resolution[2]-1) << "," <<
                    value << "\n";
            }
        }
    }

    csvfile.close();

#endif /* WITHCSVOUTPUT */
}


/* write out the grid in its actual size */
//void writeToCsvFullGrid( const Level& level ) {
void writeToCsv_full_grid( const Level& level ) {

#ifdef WITHCSVOUTPUT

    const MatrixT& grid= *level.src_grid;

    std::array< long int, 3 > corner= grid.pattern().global( {0,0,0} );

    size_t d= grid.extent(0);
    size_t h= grid.extent(1);
    size_t w= grid.extent(2);

    size_t dl= grid.local.extent(0);
    size_t hl= grid.local.extent(1);
    size_t wl= grid.local.extent(2);

    /*
    cout << "writeToCsvFullGrid " <<
        dl << "," << hl << "," << wl << " of " <<
        d << "," << h << "," << w << endl;
    */

    std::ofstream csvfile;
    std::ostringstream num_string;
    num_string << std::setw(5) << std::setfill('0') << (uint32_t) filenumber->get();
    csvfile.open( "image_unit" + std::to_string(dash::myid()) +
        ".csv." + num_string.str() );

    grid.barrier();
    if ( 0 == dash::myid() ) {
        csvfile << " z-coord,y-coord,x-coord,heat" << "\n";
        filenumber->set( 1 + (uint32_t) filenumber->get()  );
    }
    grid.barrier();

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


void initgrid( Level& level ) {

    /* not strictly necessary but it also avoids NAN values */
    dash::fill( level.src_grid->begin(), level.src_grid->end(), 5.0 );
    dash::fill( level.dst_grid->begin(), level.dst_grid->end(), 0.0 );
    dash::fill( level.rhs_grid->begin(), level.rhs_grid->end(), 0.0 );

    level.src_grid->barrier();
}


/* apply boundary value settings, where the top and bottom planes have a 
hot circle in the middle and everything else is cold */
void initboundary( Level& level ) {

    using index_t = dash::default_index_t;

    double gd= level.src_grid->extent(0);
    double gh= level.src_grid->extent(1);
    double gw= level.src_grid->extent(2);

    /* This way of setting boundaries uses subsampling on the top and bottom
    planes to determine the border values. This is another logical way that
    may be convenient sometimes. It guarantees that the boundary values on all
    the levels match.
    All other sides are constant at 0.0 degrees. The top an bottom circles are
    hot with 10.0 degrees. */

    auto lambda= [gd,gh,gw]( const auto& coords ) {

        index_t z= coords[0];
        index_t y= coords[1];
        index_t x= coords[2];

        double ret= 1.0;

        /* for simplicity make every side uniform */

        if ( -1 == z || gd == z ) {

            /* radius differs on top and bottom plane */
            //double r= ( -1 == z ) ? 0.4 : 0.3;
            double r= 0.4;
            double r2= r*r;

            double lowvalue= 2.0;
            double highvalue= 9.0;

            double midx= 0.5;
            double midy= 0.5;

            /* At entry (x/gw,y/gh) we sample the
            rectangle [ x/gw,(x+1)/gw ) x [ y/gw, (y+1)/gh ) with m² points. */
            int32_t m= 3;
            int32_t m2= m*m;

            double sum= 0.0;
            double weight= 0.0;

            for ( double iy= -m+1; iy < m; iy++ ) {
                for ( double ix= -m+1; ix < m; ix++ ) {

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

    level.src_halo->set_fixed_halos( lambda );
    level.dst_halo->set_fixed_halos( lambda );
}


/* sets all boundary values to 0, that is what is neede on the coarser grids */
void initboundary_zero( Level& level ) {

    using index_t = dash::default_index_t;

    auto lambda= []( const auto& coords ) { return 0.0; };

    level.src_halo->set_fixed_halos( lambda );
    level.dst_halo->set_fixed_halos( lambda );
}


void markunits( MatrixT& grid ) {

    /* Mark unit bordery by setting the first local rows and columns */

    size_t w= grid.local.extent(2);
    size_t h= grid.local.extent(1);
    size_t d= grid.local.extent(0);

    for ( size_t i = 0; i < d; ++i ) {
        for ( size_t j = 0; j < h; ++j ) {
            grid.local[i][j][0] = 8.0;
        }
    }

    for ( size_t i = 0; i < d; ++i ) {
        for ( size_t k = 0; k < w; ++k ) {
            grid.local[i][0][k] = 8.0;
        }
    }

    for ( size_t j = 0; j < h; ++j ) {
        for ( size_t k = 0; k < w; ++k ) {
            grid.local[0][j][k] = 8.0;
        }
    }
}


/* check some grid values for 3d mirror symmetry. This should hold for
appropriate boundary conditions and a correct solver.

Here we use global accesses for simplicity. */
bool check_symmetry( MatrixT& grid, double eps ) {

    if ( 0 == dash::myid() ) {

        size_t w= grid.extent(2);
        size_t h= grid.extent(1);
        size_t d= grid.extent(0);

        size_t m= std::min( std::min( w, h ), d ) /2;

        /* x-y-z diagonals */
        for ( size_t t= 0; t < m; ++t ) {

            double first= grid[d/2+t][h/2+t][w/2+t];

            if ( std::fabs( first - grid[d/2+t][h/2+t][w/2-t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2+t][h/2-t][w/2+t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2+t][h/2-t][w/2-t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2+t][w/2+t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2+t][w/2-t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2-t][w/2+t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2-t][w/2-t] ) > eps ) return false;
        }

        /* x-y diagonals */
        for ( size_t t= 0; t < m; ++t ) {

            double first= grid[d/2][h/2+t][w/2+t];

            if ( std::fabs( first - grid[d/2][h/2+t][w/2-t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2][h/2-t][w/2+t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2][h/2-t][w/2-t] ) > eps ) return false;
        }

        /* y-z diagonals */
        for ( size_t t= 0; t < m; ++t ) {

            double first= grid[d/2+t][h/2+t][w/2];

            if ( std::fabs( first - grid[d/2+t][h/2+t][w/2] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2+t][h/2-t][w/2] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2+t][w/2] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2-t][w/2] ) > eps ) return false;
        }

        /* x-z diagonals */
        for ( size_t t= 0; t < m; ++t ) {

            double first= grid[d/2+t][h/2][w/2+t];

            if ( std::fabs( first - grid[d/2+t][h/2][w/2-t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2][w/2+t] ) > eps ) return false;
            if ( std::fabs( first - grid[d/2-t][h/2][w/2-t] ) > eps ) return false;
        }

    }

    return true;
}


void scaledownboundary( Level& fine, Level& coarse ) {

    assert( coarse.src_grid->extent(2)*2 == fine.src_grid->extent(2) );
    assert( coarse.src_grid->extent(1)*2 == fine.src_grid->extent(1) );
    assert( coarse.src_grid->extent(0)*2 == fine.src_grid->extent(0) );

    size_t dmax= coarse.src_grid->extent(0);
    size_t hmax= coarse.src_grid->extent(1);
    //size_t wmax= coarse.src_grid->extent(2);

    auto finehalo= fine.src_halo;

    auto lambda= [&finehalo,&dmax,&hmax]( const auto& coord ) {

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

    coarse.src_halo->set_fixed_halos( lambda );
    coarse.dst_halo->set_fixed_halos( lambda );
}


void scaledown( Level& fine, Level& coarse ) {

    MatrixT& finegrid= *fine.src_grid;
    MatrixT& fine_rhs_grid= *fine.rhs_grid;
    MatrixT& coarsegrid= *coarse.src_grid;
    MatrixT& coarse_rhs_grid= *coarse.rhs_grid;

    double ax= fine.ax;
    double ay= fine.ay;
    double az= fine.az;
    double ac= fine.acenter;
    double ff= fine.ff;

    uint64_t par= finegrid.team().size();
    uint64_t param= finegrid.local.extent(0)*finegrid.local.extent(1)*finegrid.local.extent(2);

    // scaledown
    minimon.start();

    assert( (coarsegrid.extent(2)+1) * 2 == finegrid.extent(2)+1 );
    assert( (coarsegrid.extent(1)+1) * 2 == finegrid.extent(1)+1 );
    assert( (coarsegrid.extent(0)+1) * 2 == finegrid.extent(0)+1 );

    const auto& extentc= coarsegrid.pattern().local_extents();
    const auto& cornerc= coarsegrid.pattern().global( {0,0,0} );
    const auto& extentf= finegrid.pattern().local_extents();
    const auto& cornerf= finegrid.pattern().global( {0,0,0} );

    assert( cornerc[0] * 2 == cornerf[0] );
    assert( cornerc[1] * 2 == cornerf[1] );
    assert( cornerc[2] * 2 == cornerf[2] );

    assert( 0 == cornerc[0] %2 );
    assert( 0 == cornerc[1] %2 );
    assert( 0 == cornerc[2] %2 );

    assert( extentc[0] * 2 == extentf[0] || extentc[0] * 2 +1 == extentf[0] );
    assert( extentc[1] * 2 == extentf[1] || extentc[1] * 2 +1 == extentf[1] );
    assert( extentc[2] * 2 == extentf[2] || extentc[2] * 2 +1 == extentf[2] );

    /* Here we  $ r= f - Au $ on the fine grid and 'straigth injection' to the 
    rhs of the coarser grid in one. Therefore, we don't need a halo of the fine 
    grid, because the stencil neighbor points on the fine grid are always there 
    for a coarse grid point. 
    According to the text book (Introduction to Algebraic Multigrid -- Course notes 
    of an algebraic multigrid course at univertisty of Heidelberg in Wintersemester 
    1998/99, Version 1.1 by Christian Wagner http://www.mgnet.org/mgnet/papers/Wagner/amgV11.pdf) 
    there should by an extra factor 1/2^3 for the coarse value. But this doesn't seem to work, 
    factor 4.0 works much better. */

    /* this is the plain and simple and slow version. do optimization by pointer 
    arithmethic later. */
    for ( size_t z= 0; z < extentc[0] ; z++ ) {
        for ( size_t y= 0; y < extentc[1] ; y++ ) {
            for ( size_t x= 0; x < extentc[2]; x++ ) {
                coarse_rhs_grid.local[z][y][x]= 4.0 * /* extra factor? */ ( 
                    ff * fine_rhs_grid.local[2*z+1][2*y+1][2*x+1] -
                    ax * ( finegrid[2*z+1][2*y+1][2*x+0] + finegrid[2*z+1][2*y+1][2*x+2] ) -
                    ay * ( finegrid[2*z+1][2*y+0][2*x+1] + finegrid[2*z+1][2*y+2][2*x+1] ) -
                    az * ( finegrid[2*z+0][2*y+1][2*x+1] + finegrid[2*z+2][2*y+1][2*x+1] ) -
                    ac * finegrid[2*z+1][2*y+1][2*x+1]);
            }
        }
    }

    /* set coarse grid to 0.0 */
    dash::fill( coarsegrid.begin(), coarsegrid.end(), 0.0 );

    minimon.stop( "scaledown", par, param );
}


/* this version uses a correct prolongation from the coarser grid of (2^n)^3 to (2^(n+1))^3
elements. Note that it is 2^n elements per dimension instead of 2^n -1!
This version loops over the coarse grid */
//void scaleup_loop_coarse( Level& coarse, Level& fine ) {
void scaleup( Level& coarse, Level& fine ) {

    using signed_size_t = typename std::make_signed<size_t>::type;

    MatrixT& coarsegrid= *coarse.src_grid;
    MatrixT& finegrid= *fine.src_grid;

    uint64_t par= coarsegrid.team().size();
    uint64_t param= coarsegrid.local.extent(0)*coarsegrid.local.extent(1)*coarsegrid.local.extent(2);

    // scaleup
    minimon.start();

    assert( (coarsegrid.extent(2)+1) * 2 == finegrid.extent(2)+1 );
    assert( (coarsegrid.extent(1)+1) * 2 == finegrid.extent(1)+1 );
    assert( (coarsegrid.extent(0)+1) * 2 == finegrid.extent(0)+1 );

    const auto& extentc= coarsegrid.pattern().local_extents();
    const auto& cornerc= coarsegrid.pattern().global( {0,0,0} );
    const auto& extentf= finegrid.pattern().local_extents();
    const auto& cornerf= finegrid.pattern().global( {0,0,0} );

    assert( cornerc[0] * 2 == cornerf[0] );
    assert( cornerc[1] * 2 == cornerf[1] );
    assert( cornerc[2] * 2 == cornerf[2] );

    assert( 0 == cornerc[0] %2 );
    assert( 0 == cornerc[1] %2 );
    assert( 0 == cornerc[2] %2 );

    assert( extentc[0] * 2 == extentf[0] || extentc[0] * 2 +1 == extentf[0] );
    assert( extentc[1] * 2 == extentf[1] || extentc[1] * 2 +1 == extentf[1] );
    assert( extentc[2] * 2 == extentf[2] || extentc[2] * 2 +1 == extentf[2] );

    /* if last element in coarse grid per dimension has no 2*i+2 element in 
     the local fine grid, then handle it as a separate loop using halo. 
     sub[i] is always 0 or 1 */
    std::array< size_t, 3 > sub;
    for ( uint32_t i= 0; i < 3; ++i ) {
         sub[i]= ( extentc[i] * 2 == extentf[i] ) ? 1 : 0;
    }

    /* 1) start async halo exchange for coarse grid*/
    coarse.src_halo->update_async();

    /* 2) first set fine grid to 0.0, becasue afterwards there are 
    multiple += operations per fine grid element */

    /* 3) second loop over the coarse grid and add the contributions to the 
    fine grid elements */

    /* for correctness, write down the clear and simple version first and
    try to apply the optimization from scaleup_old() again -- see above.
    Make sure that this was correct though. And write down some good comment
    why it is correct! */

    for ( size_t z= 0; z < extentc[0] - sub[0]; z++ ) {
        for ( size_t y= 0; y < extentc[1] - sub[1]; y++ ) {
            for ( size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;


            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
        /* for last element when 1 == sub[1] do the same thing as above 
        but without the +2 steps in y-dimension */
        for ( size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {
            for ( size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }
    /* for last element when 1 == sub[0] do the same thing as above 
    but without the +2 steps in z-dimension */
    for ( size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {
        for ( size_t y= 0; y < extentc[1] - sub[1]; y++ ) {
            for ( size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
        /* for last element when 1 == sub[1] do the same thing as above 
        but without the +2 steps in y-dimension */
        for ( size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {
            for ( size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= coarsegrid.local[z][y][x];

                finegrid.local[2*z+1][2*y+1][2*x+1] += tmp;

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;
                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* 4) wait for async halo exchange */
    coarse.src_halo->wait();

    /* 5) do the remaining updates with contributions from the coarse halo,
    this is quite a spaghetti code for 6 sides, 12 edges, and 8 corners -- can we do this any better? */

    /* 6 planes */

    /* z= -1 */
    { signed_size_t z= -1;
        for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
        /* for last element when 1 == sub[1] do the same thing as above 
        but without the +2 steps in y-dimension */
        for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* z= N */
    if ( 0 == sub[0]) { 
        signed_size_t z= extentc[0];
        for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
        /* for last element when 1 == sub[1] do the same thing as above 
        but without the +2 steps in y-dimension */
        for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+1][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* y= -1 */
    for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {
        { signed_size_t y= -1;
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );


                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;


            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
    }
    /* for last element when 1 == sub[0] do the same thing as above 
    but without the +2 steps in z-dimension */
    for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {
        { signed_size_t y= -1;
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+2][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* y= N */
    if ( 0 == sub[1] ) {
        signed_size_t y= extentc[1];
        for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );


                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;


            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
        }

        /* for last element when 1 == sub[0] do the same thing as above 
        but without the +2 steps in z-dimension */
        for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+0][2*x+1] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* x= -1 */
    for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {
        for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;
            }
        }
        /* for last element when 1 == sub[1] do the same thing as above 
        but without the +2 steps in y-dimension */
        for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
            }
        }
    }
    /* for last element when 1 == sub[0] do the same thing as above 
    but without the +2 steps in z-dimension */
    for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {
        for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
        }
        /* for last element when 1 == sub[1] do the same thing as above 
        but without the +2 steps in y-dimension */
        for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+2] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
        }
    }

    /* x= N */
    if ( 0 == sub[2] ) {
        signed_size_t x= extentc[2];
        for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {
            for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
            /* for last element when 1 == sub[1] do the same thing as above 
            but without the +2 steps in y-dimension */
            for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
        /* for last element when 1 == sub[0] do the same thing as above 
        but without the +2 steps in z-dimension */
        for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {
            for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
            /* for last element when 1 == sub[1] do the same thing as above 
            but without the +2 steps in y-dimension */
            for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+1][2*x+0] += 0.5*tmp;

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;
                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* 3x4 edges */

    /* loop z */

    /* y= -1, x= -1 */
    { signed_size_t y= -1;
        { signed_size_t x= -1;
            for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[0] do the same thing as above 
            but without the +2 steps in z-dimension */
            for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+2][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
        }
    }

    /* y= N, x= -1 */
    if ( 0 == sub[1] ) { 
        signed_size_t y= extentc[1];
        { signed_size_t x= -1;
            for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[0] do the same thing as above 
            but without the +2 steps in z-dimension */
            for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+0][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
        }
    }

    /* y= -1, x= N */
    if ( 0 == sub[2] ) { signed_size_t x= extentc[2];
        { signed_size_t y= -1;
            for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
            /* for last element when 1 == sub[0] do the same thing as above 
            but without the +2 steps in z-dimension */
            for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+2][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* y= N, x= N */
    if ( 0 == sub[1] ) { 
        signed_size_t y= extentc[1];
        if ( 0 == sub[2] ) {
            signed_size_t x= extentc[2];
            for ( signed_size_t z= 0; z < extentc[0] - sub[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
            /* for last element when 1 == sub[0] do the same thing as above 
            but without the +2 steps in z-dimension */
            for ( signed_size_t z= extentc[0] - sub[0]; z < extentc[0]; z++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+1][2*y+0][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* loop y */

    /* z= -1, x= -1 */
    { signed_size_t z= -1;
        { signed_size_t x= -1;
            for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[1] do the same thing as above 
            but without the +2 steps in y-dimension */
            for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );


                finegrid.local[2*z+2][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
            }
        }
    }

    /* z= N, x= -1 */
    if ( 0 == sub[0] ) { 
        signed_size_t z= extentc[0];
        { signed_size_t x= -1;
            for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[1] do the same thing as above 
            but without the +2 steps in y-dimension */
            for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );


                finegrid.local[2*z+0][2*y+1][2*x+2] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
        }
    }

    /* z= -1, x= N */
    { signed_size_t z= -1;
        if ( 0 == sub[2] ) { signed_size_t x= extentc[2];
            for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
            /* for last element when 1 == sub[1] do the same thing as above 
            but without the +2 steps in y-dimension */
            for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );


                finegrid.local[2*z+2][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* z= N, x= N */
    if ( 0 == sub[0] ) { signed_size_t z= extentc[0];
        if ( 0 == sub[2] ) { signed_size_t x= extentc[2];
            for ( signed_size_t y= 0; y < extentc[1] - sub[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
            /* for last element when 1 == sub[1] do the same thing as above 
            but without the +2 steps in y-dimension */
            for ( signed_size_t y= extentc[1] - sub[1]; y < extentc[1]; y++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );


                finegrid.local[2*z+0][2*y+1][2*x+0] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }


    /* loop x */

    /* z= -1, y= -1 */
    { signed_size_t z= -1;
        { signed_size_t y= -1;
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* z= N, y= -1 */
    if ( 0 == sub[0] ) { 
        signed_size_t z= extentc[0];
        { signed_size_t y= -1;
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+2][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* z= -1, y= N */
    { signed_size_t z= -1;
        if ( 0 == sub[1] ) { signed_size_t y= extentc[1];
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+2][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* z= N, y= N */
    if ( 0 == sub[0] ) { signed_size_t z= extentc[0];
        if ( 0 == sub[1] ) { signed_size_t y= extentc[1];
            for ( signed_size_t x= 0; x < extentc[2] - sub[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
            /* for last element when 1 == sub[2] do the same thing as above 
            but without the +2 steps in x-dimension */
            for ( signed_size_t x= extentc[2] - sub[2]; x < extentc[2]; x++ ) {

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );

                finegrid.local[2*z+0][2*y+0][2*x+1] += 0.25*tmp;

                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    /* 8 corners */

    { signed_size_t z= -1;
        { signed_size_t y= -1;
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+2][2*y+2][2*x+2] += 0.125*tmp;
            }
            if ( 0 == sub[2]) { signed_size_t x= extentc[2];

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+2][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
        if ( 0 == sub[1] ) { signed_size_t y= extentc[1];
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+2][2*y+0][2*x+2] += 0.125*tmp;
            }
            if ( 0 == sub[2]) { signed_size_t x= extentc[2];

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+2][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }
    if ( 0 == sub[0] ) { signed_size_t z= extentc[0];
        { signed_size_t y= -1;
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+0][2*y+2][2*x+2] += 0.125*tmp;
            }
            if ( 0 == sub[2]) { signed_size_t x= extentc[2];

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+0][2*y+2][2*x+0] += 0.125*tmp;
            }
        }
        if ( 0 == sub[1] ) { signed_size_t y= extentc[1];
            { signed_size_t x= -1;

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+0][2*y+0][2*x+2] += 0.125*tmp;
            }
            if ( 0 == sub[2]) { signed_size_t x= extentc[2];

                double tmp= *coarse.src_halo->halo_element_at( {cornerc[0]+z,cornerc[1]+y,cornerc[2]+x} );
                finegrid.local[2*z+0][2*y+0][2*x+0] += 0.125*tmp;
            }
        }
    }

    minimon.stop( "scaleup", par, param );
}


/* this version uses a correct prolongation from the coarser grid of (2^n)^3 to (2^(n+1))^3
elements. Note that it is 2^n elements per dimension instead of 2^n -1!
This version loops over the fine grid */
void scaleup_loop_fine( Level& coarse, Level& fine ) {

    MatrixT& coarsegrid= *coarse.src_grid;
    MatrixT& finegrid= *fine.src_grid;

    uint64_t par= coarsegrid.team().size();
    uint64_t param= coarsegrid.local.extent(0)*coarsegrid.local.extent(1)*coarsegrid.local.extent(2);

    // scaleup
    minimon.start();

    assert( coarsegrid.extent(2) * 2 == finegrid.extent(2) );
    assert( coarsegrid.extent(1) * 2 == finegrid.extent(1) );
    assert( coarsegrid.extent(0) * 2 == finegrid.extent(0) );

    const auto& extentc= coarsegrid.pattern().local_extents();
    const auto& cornerc= coarsegrid.pattern().global( {0,0,0} );
    const auto& extentf= finegrid.pattern().local_extents();
    const auto& cornerf= finegrid.pattern().global( {0,0,0} );

    assert( cornerc[0] * 2 == cornerf[0] );
    assert( cornerc[1] * 2 == cornerf[1] );
    assert( cornerc[2] * 2 == cornerf[2] );

    assert( 0 == cornerc[0] %2 );
    assert( 0 == cornerc[1] %2 );
    assert( 0 == cornerc[2] %2 );

    assert( extentc[0] * 2 == extentf[0] );
    assert( extentc[1] * 2 == extentf[1] );
    assert( extentc[2] * 2 == extentf[2] );

    /* 1) start async halo exchange for coarse grid*/
    coarse.src_halo->update_async();

    /* 2) iterate over inner area of fine grid, that is leaving a border of size 2
    because the access pattern */

    for ( size_t z= 1; z < extentf[0]-1 ; z++ ) {

        //size_t dz= -1 + 2*(z&1); // same as ( 0 == z%2 ) ? -1 : +1;
        size_t dz= ( 0 == z%2 ) ? -1 : +1;

        for ( size_t y= 1; y < extentf[1]-1 ; y++ ) {
 
            //size_t dy= -1 + 2*(y&1); // same as ( 0 == y%2 ) ? -1 : +1;
            size_t dy= ( 0 == y%2 ) ? -1 : +1;

            for ( size_t x= 1; x < extentf[2]-1 ; x++ ) {

                //size_t dx= -1 + 2*(x&1); // same as ( 0 == x%2 ) ? -1 : +1;
                size_t dx= ( 0 == x%2 ) ? -1 : +1;

                finegrid.local[z][y][x]= 
                    0.75*0.75*0.75 * coarsegrid.local[z/2   ][y/2   ][x/2   ] + 
                    0.75*0.75*0.25 * coarsegrid.local[z/2   ][y/2   ][x/2+dx] + 
                    0.75*0.25*0.75 * coarsegrid.local[z/2   ][y/2+dy][x/2   ] + 
                    0.75*0.25*0.25 * coarsegrid.local[z/2   ][y/2+dy][x/2+dx] + 
                    0.25*0.75*0.75 * coarsegrid.local[z/2+dz][y/2   ][x/2   ] + 
                    0.25*0.75*0.25 * coarsegrid.local[z/2+dz][y/2   ][x/2+dx] + 
                    0.25*0.25*0.75 * coarsegrid.local[z/2+dz][y/2+dy][x/2   ] + 
                    0.25*0.25*0.25 * coarsegrid.local[z/2+dz][y/2+dy][x/2+dx];
            }
        }
    }

    /* 4) wait for async halo exchange */
    coarse.src_halo->wait();

    /* 5) do the remaining updates with contributions from the coarse halo */
    // z-y
    // y-x
    //z-x

    minimon.stop( "scaleup", par, param );
}


void transfertofewer( Level& source /* with larger team*/, Level& dest /* with smaller team */ ) {

    /* should only be called by the smaller team */
    assert( 0 == dest.src_grid->team().position() );

cout << "unit " << dash::myid() << " transfertofewer" << endl;

    /* we need to find the coordinates that the local unit needs to receive
    from several other units that are not in this team */

    /* we can safely assume that the source blocks are copied entirely */

    std::array< long int, 3 > corner= dest.src_grid->pattern().global( {0,0,0} );
    std::array< long unsigned int, 3 > sizes= dest.src_grid->pattern().local_extents();

    /*
    cout << "    start coord: " <<
        corner[0] << ", "  << corner[1] << ", " << corner[2] << endl;
    cout << "    extents: " <<
            sizes[0] << ", "  << sizes[1] << ", " << sizes[2] << endl;
    cout << "    dest local  dist " << dest.src_grid->lend() - dest.src_grid->lbegin() << endl;
    cout << "    dest global dist " << dest.src_grid->end() - dest.src_grid->begin() << endl;
    cout << "    src  local  dist " << source.src_grid->lend() - source.src_grid->lbegin() << endl;
    cout << "    src  global dist " << source.src_grid->end() - source.src_grid->begin() << endl;
    */

    /* Can I do this any cleverer than loops over the n-1 non-contiguous
    dimensions and then a dash::copy for the 1 contiguous dimension? */

    /*
    for ( uint32_t z= 0; z < sizes[0]; z++ ) {
        for ( uint32_t y= 0; y < sizes[1]; y++ ) {
            for ( uint32_t x= 0; x < sizes[2]; x++ ) {

                (*dest.src_grid)(z,y,x)= (*source.src_grid)(z,y,x);
            }
        }
    }
    */

    for ( uint32_t z= 0; z < sizes[0]; z++ ) {
        for ( uint32_t y= 0; y < sizes[1]; y++ ) {

            auto start= source.src_grid->begin() + ((corner[0]+z)*sizes[1]+y)*sizes[2];
            std::copy( start, start + sizes[2], &dest.src_grid->local[z][y][0] );

            //dash::copy( start, start + sizes[2], &dest.src_grid->local[z][y][0] );
            //dash::copy( source.grid.begin()+40, source.grid.begin()+48, buf );
        }
    }

}


void transfertomore( Level& source /* with smaller team*/, Level& dest /* with larger team */ ) {

    /* should only be called by the smaller team */
    assert( 0 == source.src_grid->team().position() );

cout << "unit " << dash::myid() << " transfertomore" << endl;

    /* we need to find the coordinates that the local unit needs to receive
    from several other units that are not in this team */

    /* we can safely assume that the source blocks are copied entirely */

    std::array< long int, 3 > corner= source.src_grid->pattern().global( {0,0,0} );
    std::array< long unsigned int, 3 > sizes= source.src_grid->pattern().local_extents();

    /*
    cout << "    start coord: " <<
        corner[0] << ", "  << corner[1] << ", " << corner[2] << endl;
    cout << "    extents: " <<
            sizes[0] << ", "  << sizes[1] << ", " << sizes[2] << endl;
    cout << "    dest local  dist " << source.src_grid->lend() - source.src_grid->lbegin() << endl;
    cout << "    dest global dist " << source.src_grid->end() - source.src_grid->begin() << endl;
    cout << "    src  local  dist " << dest.src_grid->lend() - dest.src_grid->lbegin() << endl;
    cout << "    src  global dist " << dest.src_grid->end() - dest.src_grid->begin() << endl;
    */

    /* stupid but functional version for the case with only one unit in the smaller team, very slow individual accesses */
    /*
    for ( uint32_t z= 0; z < sizes[0]; z++ ) {
        for ( uint32_t y= 0; y < sizes[1]; y++ ) {
            for ( uint32_t x= 0; x < sizes[2]; x++ ) {

                (*dest.src_grid)(z,y,x)= (*source.src_grid)(z,y,x);
            }
        }
    }
    */
    for ( uint32_t z= 0; z < sizes[0]; z++ ) {
        for ( uint32_t y= 0; y < sizes[1]; y++ ) {

            double* start= &source.src_grid->local[z][y][0];
            auto to= dest.src_grid->begin() + ((corner[0]+z)*sizes[1]+y)*sizes[2];
            std::copy( start, start + sizes[2], to );

            //dash::copy( start, start + sizes[2], &dest.src_grid->local[z][y][0] );
            //dash::copy( source.grid.begin()+40, source.grid.begin()+48, buf );
        }
    }

    //std::copy( source.src_grid->begin(), source.src_grid->end(), dest.src_grid->begin() );
}


/**
Smoothen the given level from oldgrid+src_halo to newgrid. Call Level::swap() at the end.

The parallel global residual is returned as a return parameter, but only
if it is not NULL because then the expensive parallel reduction is just avoided.
*/
double smoothen( Level& level, Allreduce& res ) {
    SCOREP_USER_FUNC()

    uint32_t par= level.src_grid->team().size();

    // smoothen
    minimon.start();

    level.src_grid->barrier();

    size_t ld= level.src_grid->local.extent(0);
    size_t lh= level.src_grid->local.extent(1);
    size_t lw= level.src_grid->local.extent(2);

    double localres= 0.0;

    double ax= level.ax;
    double ay= level.ay;
    double az= level.az;
    double ac= level.acenter;
    double ff= level.ff;
    double m= level.m;

    const double c= 1.0;

    // async halo update
    level.src_halo->update_async();

    // smoothen_inner
    minimon.start();

    // update inner

    /* the start value for both, the y loop and the x loop is 1 because either there is
    a border area next to the halo -- then the first column or row is covered below in
    the border update -- or there is an outside border -- then the first column or row
    contains the boundary values. */
    auto next_layer_off = lw * lh;
    auto core_offset = lw * (lh + 1) + 1;
    for ( size_t z= 1; z < ld-1; z++ ) {
        for ( size_t y= 1; y < lh-1; y++ ) {

            /* this should eventually be done with Alpaka or Kokkos to look
            much nicer but still be fast */

            const double* __restrict p_core=  level.src_grid->lbegin() + core_offset;
            const double* __restrict p_east=  p_core + 1;
            const double* __restrict p_west=  p_core - 1;
            const double* __restrict p_north= p_core + lw;
            const double* __restrict p_south= p_core - lw;
            const double* __restrict p_up=    p_core + next_layer_off;
            const double* __restrict p_down=  p_core - next_layer_off;
            const double* __restrict p_rhs=   level.rhs_grid->lbegin() + core_offset;
            double* __restrict p_new= level.dst_grid->lbegin() + core_offset;

            for ( size_t x= 1; x < lw-1; x++ ) {

                /*
                stability condition: r <= 1/2 with r= dt/h^2 ==> dt <= 1/2*h^2
                dtheta= ru*u_plus + ru*u_minus - 2*ru*u_center with ru=dt/hu^2 <= 1/2
                */
                double dtheta= m * (
                    ff * *p_rhs -
                    ax * ( *p_east + *p_west ) -
                    ay * ( *p_north + *p_south ) -
                    az * ( *p_up + *p_down ) -
                    ac * *p_core );
                *p_new= *p_core + c * dtheta;

                localres= std::max( localres, std::fabs( dtheta ) );

                p_core++;
                p_east++;
                p_west++;
                p_north++;
                p_south++;
                p_up++;
                p_down++;
                p_rhs++;
                p_new++;
            }
            core_offset += lw;
        }
        core_offset += 2 * lw;
    }

    minimon.stop( "smoothen_inner", par, /* elements */ (ld-2)*(lh-2)*(lw-2), /* flops */ 16*(ld-2)*(lh-2)*(lw-2), /*loads*/ 7*(ld-2)*(lh-2)*(lw-2), /* stores */ (ld-2)*(lh-2)*(lw-2) );

    // smoothen_wait
    minimon.start();
    // wait for async halo update

    level.src_halo->wait();

    minimon.stop( "smoothen_wait", par, /* elements */ ld*lh*lw );

    // smoothen_collect
    minimon.start();

    /* unit 0 (of any active team) waits until all local residuals from all
    other active units are in */
    res.collect_and_spread( level.src_grid->team() );

    minimon.stop( "smoothen_collect", par );

    // smoothen_outer
    minimon.start();

    /// begin pointer of local block, needed because halo border iterator is read-only
    auto grid_local_begin= level.dst_grid->lbegin();
    auto rhs_grid_local_begin= level.rhs_grid->lbegin();

    auto bend = level.src_halo->bend();
    // update border area
    for( auto it = level.src_halo->bbegin(); it != bend; ++it ) {

        double dtheta= m * ( 
            ff * rhs_grid_local_begin[ it.lpos() ] -
            ax * ( it.value_at(4) + it.value_at(5) ) -
            ay * ( it.value_at(2) + it.value_at(3) ) -
            az * ( it.value_at(0) + it.value_at(1) ) -
            ac * *it );
        grid_local_begin[ it.lpos() ]= *it + c * dtheta;

        localres= std::max( localres, std::fabs( dtheta ) );
    }

    minimon.stop( "smoothen_outer", par, /* elements */ 2*(ld*lh+lh*lw+lw*ld),
        /* flops */ 16*(ld*lh+lh*lw+lw*ld), /*loads*/ 7*(ld*lh+lh*lw+lw*ld), /* stores */ (ld*lh+lh*lw+lw*ld) );

    // smoothen_wait_res
    minimon.start();

    res.wait( level.src_grid->team() );

    /* global residual from former iteration */
    double oldres= res.get();

    res.set( &localres, level.src_grid->team() );

    minimon.stop( "smoothen_wait_res", par );

    level.swap();

    minimon.stop( "smoothen", par, /* elements */ ld*lh*lw,
        /* flops */ 16*ld*lh*lw, /*loads*/ 7*ld*lh*lw, /* stores */ ld*lh*lw );

    return oldres;
}

//#define DETAILOUTPUT 1

template<typename Iterator>
void v_cycle( Iterator it, Iterator itend,
        uint32_t numiter, double epsilon, Allreduce& res ) {
    SCOREP_USER_FUNC()

    Iterator itnext( it );
    ++itnext;
    /* reached end of recursion? */
    if ( itend == itnext ) {

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== start on coarsest grid ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

        /* smoothen completely  */
        uint32_t j= 0;
        res.reset( (*it)->src_grid->team() );
        while ( res.get() > epsilon ) {
            /* need global residual for iteration count */
            smoothen( **it, res );

            j++;
/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== after smoothen coarsest " << j << " times ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

        }
        if ( 0 == dash::myid()  ) {
            cout << "    smoothing coarsest " << j << " times with residual " << res.get() << endl;
        }
        writeToCsv( **it );

        return;
    }

    /* stepped on the dummy level? ... which is there to signal that it is not
    the end of the parallel recursion on the coarsest level but a subteam is
    going on to solve the coarser levels and this unit is not in that subteam.
    ... sounds complicated, is complicated, change only if you know what you
    are doing. */
    if ( NULL == *itnext ) {

        /* barrier 'Alice', belongs together with the next barrier 'Bob' below */
        (*it)->src_grid->team().barrier();

        cout << "all meet again here: I'm passive unit " << dash::myid() << endl;

        return;
    }

    /* stepped on a transfer level? */
    if ( (*it)->src_grid->team().size() != (*itnext)->src_grid->team().size() ) {

        /* only the members of the reduced team need to work, all others do siesta. */
        //if ( 0 == (*itnext)->grid.team().position() )
        assert( 0 == (*itnext)->src_grid->team().position() );
        {

            cout << "transfer to " <<
                (*it)->src_grid->extent(2) << "x" <<
                (*it)->src_grid->extent(1) << "x" <<
                (*it)->src_grid->extent(0) << " with " << (*it)->src_grid->team().size() << " units "
                " --> " <<
                (*itnext)->src_grid->extent(2) << "x" <<
                (*itnext)->src_grid->extent(1) << "x" <<
                (*itnext)->src_grid->extent(0) << " with " << (*itnext)->src_grid->team().size() << " units " << endl;

            transfertofewer( **it, **itnext );

            v_cycle( itnext, itend, numiter, epsilon, res );

            cout << "transfer back " <<
            (*itnext)->src_grid->extent(2) << "x" <<
            (*itnext)->src_grid->extent(1) << "x" <<
            (*itnext)->src_grid->extent(0) << " with " << (*itnext)->src_grid->team().size() << " units "
            " --> " <<
            (*it)->src_grid->extent(2) << "x" <<
            (*it)->src_grid->extent(1) << "x" <<
            (*it)->src_grid->extent(0) << " with " << (*it)->src_grid->team().size() << " units " <<  endl;

            transfertomore( **itnext, **it );
        }

        /* barrier 'Bob', belongs together with the previous barrier 'Alice' above */
        (*it)->src_grid->team().barrier();


        cout << "all meet again here: I'm active unit " << dash::myid() << endl;
        return;
    }


    /* **** normal recursion **** **** **** **** **** **** **** **** **** */

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== before ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

    /* smoothen fixed number of times */
    uint32_t j= 0;
    res.reset( (*it)->src_grid->team() );
    while ( res.get() > epsilon && j < numiter ) {

        /* need global residual for iteration count */
        smoothen( **it, res );

        j++;

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== after smoothen " << j << " times ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */
    }
    if ( 0 == dash::myid()  ) {
        cout << "    smoothing on way down " << j << " times with residual " << res.get() << endl;
    }

    writeToCsv( **it );

    /* scale down */
    if ( 0 == dash::myid() ) {
        cout << "scale down " <<
            (*it)->src_grid->extent(2) << "x" <<
            (*it)->src_grid->extent(1) << "x" <<
            (*it)->src_grid->extent(0) <<
            " --> " <<
            (*itnext)->src_grid->extent(2) << "x" <<
            (*itnext)->src_grid->extent(1) << "x" <<
            (*itnext)->src_grid->extent(0) << endl;
    }

    scaledown( **it, **itnext );
    writeToCsv( **itnext );

    /* recurse  */
    v_cycle( itnext, itend, numiter, epsilon, res );

    /* scale up */
    if ( 0 == dash::myid() ) {
        cout << "scale up " <<
            (*itnext)->src_grid->extent(2) << "x" <<
            (*itnext)->src_grid->extent(1) << "x" <<
            (*itnext)->src_grid->extent(0) <<
            " --> " <<
            (*it)->src_grid->extent(2) << "x" <<
            (*it)->src_grid->extent(1) << "x" <<
            (*it)->src_grid->extent(0) << endl;
    }
    scaleup( **itnext, **it );
    writeToCsv( **it );

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== back on fine grid ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */
    /* on the way up it ought to solve the grid rather completley,
    thus repeat until rey < eps here too???

    uint32_t j= 0;
    res.reset( (*it)->src_grid->team() );
    while ( res.get() > epsilon ) {

        smoothen( **it, res );

        j++;
    }
    if ( 0 == dash::myid() ) {
        cout << "smoothing on way up " << j << " times with residual " << res.get() << endl;
    }
    */

    j= 0;
    res.reset( (*it)->src_grid->team() );
    while ( res.get() > epsilon && j < numiter ) {

        /* need global residual for iteration count */
        smoothen( **it, res );

        j++;

        if ( 0 == dash::myid() && 0 == j % 1000 ) {

            cout << "        j " << j << " res " << res.get() << endl;
        }

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== final smoothen " << j << " times ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

    }
    if ( 0 == dash::myid() ) {
        cout << "    smoothing on way up " << j << " times with residual " << res.get() << endl;
    }

    writeToCsv( **it );
}


template<typename Iterator>
void w_cycle( Iterator it, Iterator itend,
        uint32_t numiter, double epsilon, Allreduce& res ) {
    SCOREP_USER_FUNC()

    Iterator itnext( it );
    ++itnext;
    /* reached end of recursion? */
    if ( itend == itnext ) {

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== start on coarsest grid ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

        /* smoothen completely  */
        uint32_t j= 0;
        res.reset( (*it)->src_grid->team() );
        while ( res.get() > epsilon ) {
            /* need global residual for iteration count */
            smoothen( **it, res );

            j++;
/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== after smoothen coarsest " << j << " times ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

        }
        if ( 0 == dash::myid()  ) {
            cout << "    smoothing coarsest " << j << " times with residual " << res.get() << endl;
        }
        writeToCsv( **it );

        return;
    }

    /* stepped on the dummy level? ... which is there to signal that it is not
    the end of the parallel recursion on the coarsest level but a subteam is
    going on to solve the coarser levels and this unit is not in that subteam.
    ... sounds complicated, is complicated, change only if you know what you
    are doing. */
    if ( NULL == *itnext ) {

        /* barrier 'Alice', belongs together with the next barrier 'Bob' below */
        (*it)->src_grid->team().barrier();

        cout << "all meet again here: I'm passive unit " << dash::myid() << endl;

        return;
    }

    /* stepped on a transfer level? */
    if ( (*it)->src_grid->team().size() != (*itnext)->src_grid->team().size() ) {

        /* only the members of the reduced team need to work, all others do siesta. */
        //if ( 0 == (*itnext)->grid.team().position() )
        assert( 0 == (*itnext)->src_grid->team().position() );
        {

            cout << "transfer to " <<
                (*it)->src_grid->extent(2) << "x" <<
                (*it)->src_grid->extent(1) << "x" <<
                (*it)->src_grid->extent(0) << " with " << (*it)->src_grid->team().size() << " units "
                " --> " <<
                (*itnext)->src_grid->extent(2) << "x" <<
                (*itnext)->src_grid->extent(1) << "x" <<
                (*itnext)->src_grid->extent(0) << " with " << (*itnext)->src_grid->team().size() << " units " << endl;

            transfertofewer( **it, **itnext );

            w_cycle( itnext, itend, numiter, epsilon, res );

            cout << "transfer back " <<
            (*itnext)->src_grid->extent(2) << "x" <<
            (*itnext)->src_grid->extent(1) << "x" <<
            (*itnext)->src_grid->extent(0) << " with " << (*itnext)->src_grid->team().size() << " units "
            " --> " <<
            (*it)->src_grid->extent(2) << "x" <<
            (*it)->src_grid->extent(1) << "x" <<
            (*it)->src_grid->extent(0) << " with " << (*it)->src_grid->team().size() << " units " <<  endl;

            transfertomore( **itnext, **it );
        }

        /* barrier 'Bob', belongs together with the previous barrier 'Alice' above */
        (*it)->src_grid->team().barrier();


        cout << "all meet again here: I'm active unit " << dash::myid() << endl;
        return;
    }


    /* **** normal recursion **** **** **** **** **** **** **** **** **** */

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== before ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

    /* smoothen fixed number of times */
    uint32_t j= 0;
    res.reset( (*it)->src_grid->team() );
    while ( res.get() > epsilon && j < numiter ) {

        /* need global residual for iteration count */
        smoothen( **it, res );

        j++;

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== after smoothen " << j << " times ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */
    }
    if ( 0 == dash::myid()  ) {
        cout << "    smoothing on way down " << j << " times with residual " << res.get() << endl;
    }

    writeToCsv( **it );

    /* scale down */
    if ( 0 == dash::myid() ) {
        cout << "scale down " <<
            (*it)->src_grid->extent(2) << "x" <<
            (*it)->src_grid->extent(1) << "x" <<
            (*it)->src_grid->extent(0) <<
            " --> " <<
            (*itnext)->src_grid->extent(2) << "x" <<
            (*itnext)->src_grid->extent(1) << "x" <<
            (*itnext)->src_grid->extent(0) << endl;
    }

    scaledown( **it, **itnext );
    writeToCsv( **itnext );

    /* recurse first */
    w_cycle( itnext, itend, numiter, epsilon, res );

    j= 0;
    res.reset( (*it)->src_grid->team() );
    while ( res.get() > epsilon && j < numiter ) {

        /* need global residual for iteration count */
        smoothen( **it, res );

        j++;

        if ( 0 == dash::myid() && 0 == j % 1000 ) {

            cout << "        j " << j << " res " << res.get() << endl;
        }

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== final smoothen " << j << " times ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

    }
    if ( 0 == dash::myid() ) {
        cout << "    smoothing in the middle " << j << " times with residual " << res.get() << endl;
    }


    /* recurse second */
    w_cycle( itnext, itend, numiter, epsilon, res );

    /* scale up */
    if ( 0 == dash::myid() ) {
        cout << "scale up " <<
            (*itnext)->src_grid->extent(2) << "x" <<
            (*itnext)->src_grid->extent(1) << "x" <<
            (*itnext)->src_grid->extent(0) <<
            " --> " <<
            (*it)->src_grid->extent(2) << "x" <<
            (*it)->src_grid->extent(1) << "x" <<
            (*it)->src_grid->extent(0) << endl;
    }
    scaleup( **itnext, **it );
    writeToCsv( **it );

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== back on fine grid ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */
    /* on the way up it ought to solve the grid rather completley,
    thus repeat until rey < eps here too???

    uint32_t j= 0;
    res.reset( (*it)->src_grid->team() );
    while ( res.get() > epsilon ) {

        smoothen( **it, res );

        j++;
    }
    if ( 0 == dash::myid() ) {
        cout << "smoothing on way up " << j << " times with residual " << res.get() << endl;
    }
    */

    j= 0;
    res.reset( (*it)->src_grid->team() );
    while ( res.get() > epsilon && j < numiter ) {

        /* need global residual for iteration count */
        smoothen( **it, res );

        j++;

        if ( 0 == dash::myid() && 0 == j % 1000 ) {

            cout << "        j " << j << " res " << res.get() << endl;
        }

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== final smoothen " << j << " times ==" << endl;
    cout << "  src_grid" << endl;
    (*it)->printout();
    cout << "  rhs_grid" << endl;
    (*it)->printout_rhs();
}
#endif /* DETAILOUTPUT */

    }
    if ( 0 == dash::myid() ) {
        cout << "    smoothing on way up " << j << " times with residual " << res.get() << endl;
    }

    writeToCsv( **it );
}


void smoothen_final( Level& level, double epsilon, Allreduce& res ) {
    SCOREP_USER_FUNC()

    uint64_t par= level.src_grid->team().size() ;

    // smooth_final
    minimon.start();

    uint32_t j= 0;
    res.reset( level.src_grid->team() );
    while ( res.get() > epsilon ) {

        smoothen( level, res );
        j++;
        if ( ( 0 == dash::myid() ) && ( 0 == j % 100 ) ) {
            cout << j << " smoothen finest, residual " << res.get() << "    "<< fflush << "\r";
        }

    }
    if ( 0 == dash::myid() ) {
        cout << "smoothing: " << j << " steps finest with residual " << res.get() << endl;
    }

/////////////////////////////////////////
#ifdef DETAILOUTPUT
if ( 0 == dash::myid()  ) {
    cout << "== after final smoothing ==" << endl;
    cout << "  src_grid" << endl;
    level.printout();
    cout << "  rhs_grid" << endl;
    level.printout_rhs();
}
#endif /* DETAILOUTPUT */

    minimon.stop( "smooth_final", par );
}


void do_multigrid_iteration( uint32_t howmanylevels, double eps ) {
    SCOREP_USER_FUNC()

    // setup
    minimon.start();

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    /* determine factors for width and height such that every unit has a power of two
    extent in every dimension and that the area is close to a square
    with aspect ratio \in [0,75,1,5] */

    vector<Level*> levels;
    levels.reserve( howmanylevels );

    if ( 0 == dash::myid() ) {

        cout << "run multigrid iteration with " << dash::Team::All().size() << " units "
            "for with grids from " <<
            2 << "x" <<
            2 << "x" <<
            2 <<
            " to " <<
            ((1<<(howmanylevels))-1) << "x" <<
            ((1<<(howmanylevels))-1) << "x" <<
            ((1<<(howmanylevels))-1) <<
            endl;
    }

    /* finest grid needs to be larger than 2*teamspec per dimension, 
    that means local grid is >= 2 elements */
    assert( (1<<(howmanylevels))-1 >= 2*teamspec.num_units(0) );
    assert( (1<<(howmanylevels))-1 >= 2*teamspec.num_units(1) );
    assert( (1<<(howmanylevels))-1 >= 2*teamspec.num_units(2) );

    /* create all grid levels, starting with the finest and ending with 2x2,
    The finest level is outside the loop because it is always done by dash::Team::All() */

    if ( 0 == dash::myid() ) {
        cout << "finest level is " <<
            (1<<(howmanylevels))-1 << "x" <<
            (1<<(howmanylevels))-1 << "x" <<
            (1<<(howmanylevels))-1 <<
            " distributed over " <<
            teamspec.num_units(0) << "x" <<
            teamspec.num_units(1) << "x" <<
            teamspec.num_units(2) << " units" << endl;
    }

    levels.push_back( new Level( 1.0, 1.0, 1.0,
        (1<<(howmanylevels))-1,
        (1<<(howmanylevels))-1,
        (1<<(howmanylevels))-1,
        dash::Team::All(), teamspec ) );

    /* only do initgrid on the finest level, use scaledownboundary for all others */
    initboundary( *levels.back() );
    sanitycheck( *levels.back()->src_grid );

    dash::barrier();

    --howmanylevels;
    while ( (1<<(howmanylevels))-1 >= 2*teamspec.num_units(0) && 
            (1<<(howmanylevels))-1 >= 2*teamspec.num_units(1) &&
            (1<<(howmanylevels))-1 >= 2*teamspec.num_units(2) ) {

        /*
        if ( 0 == dash::myid() ) {
            cout << "compute level " << l << " is " <<
                (1<<(howmanylevels))-1 << "x" <<
                (1<<(howmanylevels))-1 << "x" <<
                (1<<(howmanylevels))-1 <<
                " distributed over " <<
                teamspec.num_units(0) << "x" <<
                teamspec.num_units(1) << "x" <<
                teamspec.num_units(2) << " units" << endl;
        }
        */

        /* do not try to allocate >= 8GB per core -- try to prevent myself
        from running too big a simulation on my laptop */
        assert( ((1<<(howmanylevels))-1) *
            ((1<<(howmanylevels))-1) *
            ((1<<(howmanylevels))-1) < dash::Team::All().size() * (1<<27) );

        Level& previouslevel= *levels.back();

        levels.push_back(
            new Level( previouslevel,
                       (1<<(howmanylevels))-1,
                       (1<<(howmanylevels))-1,
                       (1<<(howmanylevels))-1,
                       dash::Team::All(), teamspec ) );

        /* scaledown boundary instead of initializing it from the same
        procedure, because this is very prone to subtle mistakes which
        makes the entire multigrid algorithm misbehave. */
        //scaledownboundary( previouslevel, *levels.back() );

        initboundary_zero( *levels.back() );

        dash::barrier();
        --howmanylevels;
    }

    /* here all units and all teams meet again, those that were active for the coarsest
    levels and those that were dormant */
    dash::Team::All().barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here
    but we do it for demonstration in the graphical output */
    initgrid( *levels.front() );
    //markunits( *levels.front()->src_grid );

    writeToCsv( *levels.front() );

    dash::Team::All().barrier();

    Allreduce res( dash::Team::All() );

    minimon.stop( "setup", dash::Team::All().size() );
/*
    if ( 0 == dash::myid()  ) {
        cout << endl << "start v-cycle with res " << 0.1 << endl << endl;
    }
    v_cycle( levels.begin(), levels.end(), 2, 0.1, res );
    dash::Team::All().barrier();

*/
    int n= 20;
for ( int v=0; v < 2; ++v ) {
    if ( 0 == dash::myid()  ) {
        cout << endl << "start v-cycle with res " << eps << endl << endl;
    }
    v_cycle( levels.begin(), levels.end(), n, eps, res );
    dash::Team::All().barrier();
}
/*
    if ( 0 == dash::myid()  ) {
        cout << endl << "start w-cycle with res " << eps << endl << endl;
    }
    w_cycle( levels.begin(), levels.end(), n, eps, res );
    dash::Team::All().barrier();
*/

    if ( 0 == dash::myid()  ) {
        cout << endl << "final smoothing with res " << eps << endl << endl;
    }
    smoothen_final( *levels.front(), eps, res );
    writeToCsv( *levels.front() );

    dash::Team::All().barrier();

    if ( 0 == dash::myid() ) {

        if ( ! check_symmetry( *levels.front()->src_grid, 0.001 ) ) {

            cout << "test for asymmetry of soution failed!" << endl;
        }
    }
}


/* elastic mode runs but still seems to have errors in it */
void do_multigrid_elastic( uint32_t howmanylevels, double eps ) {

    // setup
    minimon.start();

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    /* determine factors for width and height such that every unit has a power of two
    extent in every dimension and that the area is close to a square
    with aspect ratio \in [0,75,1,5] */

    uint32_t factor_z= 1;
    uint32_t factor_y= 1;
    uint32_t factor_x= 1;

    vector<Level*> levels;
    levels.reserve( howmanylevels );

    if ( 0 == dash::myid() ) {

        cout << "run elastic multigrid iteration with " << dash::Team::All().size() << " units "
            "for with grids from " <<
            2*factor_z << "x" <<
            2*factor_y << "x" <<
            2* factor_x <<
            " to " <<
            ((1<<(howmanylevels))-1)*factor_z << "x" <<
            ((1<<(howmanylevels))-1)*factor_y << "x" <<
            ((1<<(howmanylevels))-1)*factor_x <<
            endl << endl;
    }

    /* create all grid levels, starting with the finest and ending with 2x2,
    The finest level is outside the loop because it is always done by dash::Team::All() */

    if ( 0 == dash::myid() ) {
        cout << "finest level is " <<
            ((1<<(howmanylevels))-1)*factor_z << "x" <<
            ((1<<(howmanylevels))-1)*factor_y << "x" <<
            ((1<<(howmanylevels))-1)*factor_x <<
            " distributed over " <<
            teamspec.num_units(0) << "x" <<
            teamspec.num_units(1) << "x" <<
            teamspec.num_units(2) << " units" << endl;
    }

    levels.push_back( new Level( 1.0, 1.0, 1.0, 
        ((1<<(howmanylevels))-1)*factor_z ,
        ((1<<(howmanylevels))-1)*factor_y ,
        ((1<<(howmanylevels))-1)*factor_x ,
        dash::Team::All(), teamspec ) );

    /* only do initgrid on the finest level, use caledownboundary for all others */
    initboundary( *levels.back() );
    sanitycheck( *levels.back()->src_grid );

    dash::barrier();

    --howmanylevels;
    while ( 0 < howmanylevels ) {

        dash::Team& previousteam= levels.back()->src_grid->team();
        dash::Team& currentteam= ( 4 == howmanylevels ) ? previousteam.split(2) : previousteam;
        TeamSpecT localteamspec( currentteam.size(), 1, 1 );
        localteamspec.balance_extents();

        /* this is the real iteration condition for this loop! */
        if ( (1<<(howmanylevels))-1 < 2*localteamspec.num_units(0) ||
                (1<<(howmanylevels))-1 < 2*localteamspec.num_units(1) ||
                (1<<(howmanylevels))-1 < 2*localteamspec.num_units(2) ) break;

        if ( 0 == currentteam.position() ) {

            if ( previousteam.size() != currentteam.size() ) {

                /* the team working on the following grid layers has just
                been reduced. Therefore, we add an additional grid with the
                same size as the previous one but for the reduced team. Then,
                copying the data from the domain of the larger team to the
                domain of the smaller team is easy. */

                /*
                if ( 0 == currentteam.myid() ) {
                    cout << "transfer level " << 
                        ((1<<(howmanylevels+1))-1)*factor_z << "x" <<
                        ((1<<(howmanylevels+1))-1)*factor_y << "x" <<
                        ((1<<(howmanylevels+1))-1)*factor_x <<
                        " distributed over " <<
                        localteamspec.num_units(0) << "x" <<
                        localteamspec.num_units(1) << "x" <<
                        localteamspec.num_units(2) << " units" << endl;
                }
                */
                levels.push_back(
                    new Level( 1.0, 1.0, 1.0,
                               ((1<<(howmanylevels+1))-1)*factor_z,
                               ((1<<(howmanylevels+1))-1)*factor_y,
                               ((1<<(howmanylevels+1))-1)*factor_x,
                               currentteam, localteamspec ) );
                initboundary_zero( *levels.back() );
            }

            //cout << "working unit " << dash::myid() << " / " << currentteam.myid() << " in subteam at position " << currentteam.position() << endl;

            /*
            if ( 0 == currentteam.myid() ) {
                cout << "compute level " << 
                    ((1<<(howmanylevels))-1)*factor_z << "x" <<
                    ((1<<(howmanylevels))-1)*factor_y << "x" <<
                    ((1<<(howmanylevels))-1)*factor_x <<
                    " distributed over " <<
                    localteamspec.num_units(0) << "x" <<
                    localteamspec.num_units(1) << "x" <<
                    localteamspec.num_units(2) << " units" << endl;
            }
            */

            /* do not try to allocate >= 8GB per core -- try to prevent myself
            from running too big a simulation on my laptop */
            assert( ((1<<(howmanylevels))-1)*factor_z *
                    ((1<<(howmanylevels))-1)*factor_y *
                    ((1<<(howmanylevels))-1)*factor_x < currentteam.size() * (1<<27) );

            levels.push_back(
                new Level( 1.0, 1.0, 1.0,
                           ((1<<(howmanylevels))-1)*factor_z ,
                           ((1<<(howmanylevels))-1)*factor_y ,
                           ((1<<(howmanylevels))-1)*factor_x ,
                           currentteam, localteamspec ) );

            initboundary_zero( *levels.back() );

        } else {

            //cout << "waiting unit " << dash::myid() << " / " << currentteam.myid() << " in subteam at position " << currentteam.position() << endl;

            /* this is a passive unit not takin part in the subteam that
            handles the coarser grids. insert a dummy entry in the vector
            of levels to signal that this is not the coarsest level globally. */
            levels.push_back( NULL );

            break;
        }

        --howmanylevels;
    }

    /* here all units and all teams meet again, those that were active for the coarsest
    levels and those that were dormant */
    dash::Team::All().barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here
    but we do it for demonstration in the graphical output */
    initgrid( *levels.front() );
    markunits( *levels.front()->src_grid );

    writeToCsv( *levels.front() );

    dash::Team::All().barrier();

    Allreduce res( dash::Team::All() );

    minimon.stop( "setup", dash::Team::All().size() );
/*
    if ( 0 == dash::myid()  ) {
        cout << endl << "start v-cycle with res " << 0.1 << endl << endl;
    }
    v_cycle( levels.begin(), levels.end(), 2, 0.1, res );
    dash::Team::All().barrier();


    if ( 0 == dash::myid()  ) {
        cout << endl << "start v-cycle with res " << 0.01 << endl << endl;
    }
    v_cycle( levels.begin(), levels.end(), 2, 0.01, res );
    dash::Team::All().barrier();
*/

    if ( 0 == dash::myid()  ) {
        cout << endl << "start v-cycle with res " << eps << endl << endl;
    }
    v_cycle( levels.begin(), levels.end(), 4, eps, res );
    dash::Team::All().barrier();


    writeToCsv( *levels.front() );

    dash::Team::All().barrier();

    if ( 0 == dash::myid() ) {

        if ( ! check_symmetry( *levels.front()->src_grid, 0.001 ) ) {

            cout << "test for asymmetry of soution failed!" << endl;
        }
    }
}


void do_flat_iteration( uint32_t howmanylevels, double eps ) {


    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    /* determine factors for width and height such that every unit has a power of two
    extent in every dimension and that the area is close to a square
    with aspect ratio \in [0,75,1,5] */

    uint32_t factor_z= 1;
    uint32_t factor_y= 1;
    uint32_t factor_x= 1;

    if ( 0 == dash::myid() ) {

        cout << "run flat iteration with " << dash::Team::All().size() << " units "
            "for grids of " <<
            ((1<<(howmanylevels))-1)*factor_z << "x" <<
            ((1<<(howmanylevels))-1)*factor_y << "x" <<
            ((1<<(howmanylevels))-1)*factor_x <<
            endl << endl;
    }

    Level* level= new Level( 1.0, 1.0, 1.0,
        ((1<<(howmanylevels))-1)*factor_z ,
        ((1<<(howmanylevels))-1)*factor_y ,
        ((1<<(howmanylevels))-1)*factor_x ,
        dash::Team::All(), teamspec );

    dash::barrier();

    initboundary( *level );

    initgrid( *level );
    //markunits( *level->src_grid );

    writeToCsv( *level );

    dash::barrier();


    // flat_iteration
    minimon.start();

    Allreduce res( dash::Team::All() );

    uint32_t j= 0;
    while ( res.get() > eps && j < 100000 ) {

        smoothen( *level, res );

        if ( 0 == j % 100 ) {
            writeToCsv( *level );
        }

        j++;

        if ( 0 == dash::myid() && ( 0 == j % 100 ) ) {
            cout << j << ": smoothen grid with residual " << res.get() << endl;
        }
    }
    writeToCsv( *level );

    minimon.stop( "flat_iteration", dash::Team::All().size() );

    if ( 0 == dash::myid() ) {

        if ( ! check_symmetry( *level->src_grid, 0.01 ) ) {

            cout << "test for asymmetry of soution failed!" << endl;
        }
    }

    delete level;
    level= NULL;
}


bool do_test_old_scaledown() {

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    Level* a= new Level( 1.0, 1.0, 1.0, 15, 15, 15, dash::Team::All(), teamspec );
    Level* b= new Level( 1.0, 1.0, 1.0, 7, 7, 7, dash::Team::All(), teamspec );

    dash::fill( a->src_grid->begin(), a->src_grid->end(), 1 );
    a->src_grid->barrier();
    dash::fill( b->src_grid->begin(), b->src_grid->end(), 0 );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        a->printout();
        //b->printout();
    }

    b->src_grid->barrier();
    scaledown( *a, *b );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        b->printout();
    }

    b->src_grid->barrier();

    delete a;
    delete b;

    return true;
}


bool do_test_old_scaleup() {

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    Level* a= new Level( 1.0, 1.0, 1.0, 4, 4, 4, dash::Team::All(), teamspec );
    Level* b= new Level( 1.0, 1.0, 1.0, 8, 8, 8, dash::Team::All(), teamspec );

    dash::fill( a->src_grid->begin(), a->src_grid->end(), 1 );
    a->src_grid->barrier();
    dash::fill( b->src_grid->begin(), b->src_grid->end(), 0 );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        a->printout();
        //b->printout();
    }

    b->src_grid->barrier();
    scaleup( *a, *b );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        b->printout();
    }

    b->src_grid->barrier();

    delete a;
    delete b;

    return true;
}


bool do_test_new_scaleup_loop_coarse() {

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    Level* a= new Level( 1.0, 1.0, 1.0, 7, 7, 7, dash::Team::All(), teamspec );
    Level* b= new Level( 1.0, 1.0, 1.0, 15, 15, 15, dash::Team::All(), teamspec );

    dash::fill( a->src_grid->begin(), a->src_grid->end(), 1 );
    a->src_grid->barrier();
    a->src_halo->set_fixed_halos(  []( const auto& coords ) { return 1.0; } );
    dash::fill( b->src_grid->begin(), b->src_grid->end(), 0 );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        a->printout();
        //b->printout();
    }

    b->src_grid->barrier();
    //scaleup_loop_coarse( *a, *b );
    scaleup( *a, *b );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        b->printout();
    }

    b->src_grid->barrier();

    delete a;
    delete b;

    return true;
}


bool do_test_new_scaleup_loop_fine() {

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    Level* a= new Level( 1.0, 1.0, 1.0, 7, 7, 7, dash::Team::All(), teamspec );
    Level* b= new Level( 1.0, 1.0, 1.0, 15, 15, 15, dash::Team::All(), teamspec );

    dash::fill( a->src_grid->begin(), a->src_grid->end(), 1.0 );
    a->src_grid->barrier();
    dash::fill( b->src_grid->begin(), b->src_grid->end(), 0.0 );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        a->printout();
        //b->printout();
    }

    dash::barrier();
    scaleup_loop_fine( *a, *b );
    dash::barrier();

    if ( 0 == dash::myid() ) {

        b->printout();
    }

    b->src_grid->barrier();

    delete a;
    delete b;

    return true;
}


bool do_test_initboundary() {

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    //Level* a= new Level( 1.0, 1.0, 1.0, 15, 15, 15, dash::Team::All(), teamspec );
    Level* a= new Level( 1.0, 1.0, 1.0, 7, 7, 7, dash::Team::All(), teamspec );

    initboundary( *a );

    a->printout_halo();

    delete a;

    return true;
}

bool do_test_writetocsv() {

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    Level* a= new Level( 1.0, 1.0, 1.0, 15, 15, 15, dash::Team::All(), teamspec );

    initboundary( *a );
    initgrid( *a );

    writeToCsv( *a );

    delete a;

    return true;
}


bool do_test_scaleupdown() {

    if ( 0 == dash::myid() ) {
        cout << "== test scaleupdown ==" << endl;
    }

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1 );
    teamspec.balance_extents();

    Level* a= new Level( 1.0, 1.0, 1.0, 15, 15, 15, dash::Team::All(), teamspec );
    Level* b= new Level( 1.0, 1.0, 1.0, 7, 7, 7, dash::Team::All(), teamspec );

    /* fill scr_grid and _dst_grid such that it looks like the 
    residual is 0.1 in every element from a previous smoothen step */
    dash::fill( a->src_grid->begin(), a->src_grid->end(), 1 );
    dash::fill( a->dst_grid->begin(), a->dst_grid->end(), 0.9 );

    if ( 0 == dash::myid() ) {
        cout << "fill a->src_grid with 1.0 and a->dst_grid with 0.9" << endl;
    }

    a->src_grid->barrier();

    /* scale down, b should be 0.1 everywhere */
    scaledown( *a, *b );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        cout << "Scale down a to b, should produce 0.1 in b->rhs" << endl;
        b->printout_rhs();
    }

    b->src_grid->barrier();

    Allreduce res( dash::Team::All() );
    res.reset( dash::Team::All() );
    smoothen( *b, res );
    b->src_grid->barrier();

    if ( 0 == dash::myid() ) {

        cout << "smoothen once on coarse grid" << endl;
        b->printout();
    }

    scaleup( *b, *a );
    a->src_grid->barrier();


    if ( 0 == dash::myid() ) {

        cout << "Scale up adding to a, should produce 1.1 in the inner area of a" << endl;
        a->printout();
    }

    a->src_grid->barrier();

    delete a;
    delete b;

    return true;
}


/* a number of tests ... restructure the code to headers and source files, then 
do a separate main() with all the tests in a separate source file */
bool do_tests( ) {

    bool success;
    bool allsuccess= true;

    if ( 0 == dash::myid() ) { cout << "run built-in tests:" << endl; }

    /*
    success= do_test_old_scaledown();
    if ( 0 == dash::myid() ) { cout << "   old scaledown: " << success << endl; }
    allsuccess &= success;
    */

    /*
    success= do_test_old_scaleup();
    if ( 0 == dash::myid() ) { cout << "   old scaleup: " << success << endl; }
    allsuccess &= success;
    */

    /*
    success= do_test_new_scaleup_loop_coarse();
    if ( 0 == dash::myid() ) { cout << "   new scaleup loop coarse: " << success << endl; }
    allsuccess &= success;
    */

    /*
    success= do_test_new_scaleup_loop_fine();
    if ( 0 == dash::myid() ) { cout << "   new scaleup loop coarse: " << success << endl; }
    allsuccess &= success;
    */

    success= do_test_initboundary();
    if ( 0 == dash::myid() ) { cout << "   initboundary: " << success << endl; }
    allsuccess &= success;

    //success= do_test_writetocsv();
    //if ( 0 == dash::myid() ) { cout << "   writetocsv: " << success << endl; }
    //allsuccess &= success;

    /*
     * success= do_test_scaleupdown();
    if ( 0 == dash::myid() ) { cout << "   new scaleupdown: " << success << endl; }
    allsuccess &= success;
    */

    if ( 0 == dash::myid() ) { cout << "all tests: " << allsuccess << endl; }

    return allsuccess;
}


int main( int argc, char* argv[] ) {

    // main
    minimon.start();

    // dash::init
    minimon.start();
    dash::init(&argc, &argv);
    auto id= dash::myid();
    minimon.stop( "dash::init", dash::Team::All().size() );

#ifdef WITHCSVOUTPUT
    filenumber= new dash::Shared<uint32_t>();
    if ( 0 == dash::myid() ) {
        filenumber->set( 0 );
    }
    filenumber->barrier();
#endif /* WITHCSVOUTPUT */

    bool do_alltests= false;
    bool do_flatsolver= false;
    bool do_elastic= false;
    uint32_t howmanylevels= 5;
    double epsilon= 1.0e-3;

    for ( int a= 1; a < argc; a++ ) {

        if ( 0 == strncmp( "-h", argv[a], 2  ) ||
                0 == strncmp( "--help", argv[a], 6 ) ) {

            if ( 0 == dash::myid() ) {

                cout << "call me as [mpirun] '" << argv[0] << "' [-h|--help] [number-of-levels=5]" << endl;
            }
            exit(0);

        } else if ( 0 == strncmp( "-t", argv[a], 2  ) ||
                0 == strncmp( "--tests", argv[a], 7 )) {

            do_alltests= true;
            if ( 0 == dash::myid() ) {

                cout << "run tests" << endl;
            }

        } else if ( 0 == strncmp( "-f", argv[a], 2  ) ||
                0 == strncmp( "--flat", argv[a], 6 )) {

            do_flatsolver= true;
            if ( 0 == dash::myid() ) {

                cout << "do flat iteration instead of multigrid" << endl;
            }

        } else if ( 0 == strncmp( "-e", argv[a], 2  ) ||
                0 == strncmp( "--elastic", argv[a], 9 )) {

            do_elastic= true;
            if ( 0 == dash::myid() ) {

                cout << "do multigrid iteration with changing number of units per grid" << endl;
            }

        } else if ( 0 == strncmp( "--eps", argv[a], 5  ) && ( a+1 < argc ) ) {

            epsilon= atof( argv[a+1] );
            ++a;
            if ( 0 == dash::myid() ) {

                cout << "using epsilon " << epsilon << endl;
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

    if ( do_alltests ) {

        do_tests();

    } else if ( do_flatsolver ) {

        do_flat_iteration( howmanylevels, epsilon );

    } else if ( do_elastic ) {

        do_multigrid_elastic( howmanylevels, epsilon );

    } else {

        do_multigrid_iteration( howmanylevels, epsilon );
    }

#ifdef WITHCSVOUTPUT

    delete filenumber;

#endif /* WITHCSVOUTPUT */

    // dash::finalize
    minimon.start();

    dash::finalize();

    minimon.stop( "dash::finalize", dash::Team::All().size() );

    minimon.stop( "main", dash::Team::All().size() );
    minimon.print(id);

    return 0;
}
