#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>

#include <libdash.h>
#include <dash/experimental/HaloMatrixWrapper.h>

#define WITHCSVOUTPUT 1
#ifdef WITHCSVOUTPUT

uint32_t filenumber= 0;
const std::string filename("image.csv");

#endif /* WITHCSVOUTPUT */

using std::cout;
using std::endl;
using std::setw;
using std::vector;

using TeamSpecT = dash::TeamSpec<3>;
using MatrixT = dash::NArray<double,3>;
using StencilT = dash::experimental::Stencil<3>;
using StencilSpecT = dash::experimental::StencilSpec<3,6>;
using HaloMatrixWrapperT = dash::experimental::HaloMatrixWrapper<MatrixT, StencilSpecT>;

const StencilSpecT stencil_spec
({
    StencilT(-1,0,0), StencilT(1,0,0),
    StencilT(0,-1,0), StencilT(0,1,0),
    StencilT(0,0,-1), StencilT(0,0,1)
});

struct Level
{
    MatrixT grid;
    HaloMatrixWrapperT halo;

    Level( size_t w, size_t h, size_t d, dash::Team& team, TeamSpecT teamspec )
    : grid( dash::SizeSpec<3>( h+2, w+2, d+2),
            dash::DistributionSpec<3>( dash::BLOCKED, dash::BLOCKED, dash::BLOCKED ), team, teamspec ),
      halo( grid, stencil_spec )
      {}
      Level() = delete;

      ~Level() {}
};

#ifdef WITHCSVOUTPUT
void toCSV( const MatrixT& cube )
{
    if ( 0 != dash::myid() ) return;

    size_t d= cube.extent(2);
    size_t w= cube.extent(1);
    size_t h= cube.extent(0);

    std::ofstream csvfile;
    csvfile.open( filename + '.' + std::to_string(filenumber++) );
    csvfile << "x coord, y coord, z coord, heat" << "\n";
    for ( size_t z = 0; z < d; ++z )
    {
        for ( size_t y = 0; y < h; ++y )
        {
            for ( size_t x = 0; x < w; ++x )
            {
                csvfile << x << "," << y << "," << z << "," << (double) cube[y][x][z] << "\n";
            }
        }
    }
    csvfile.close();
}
#endif /* WITHCSVOUTPUT */

void writeToCsv( const MatrixT& cube )
{
#ifdef WITHCSVOUTPUT
    cube.barrier();

    if(0 == dash::myid())
        toCSV( cube );

    cube.barrier();
#endif /* WITHCSVOUTPUT */
}

void sanitycheck( const MatrixT& grid )
{
    /* check if the sum of the local extents of the matrix blocks sum up to
     * the global extents, abort otherwise */
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
                grid.extent(0) * grid.extent(1) * grid.extent(2) <<
                " == " << grid.extent(0) << "x" << grid.extent(1) << "x" << grid.extent(2) <<
                " but sum of local sizes is " << (size_t) sums[0] << std::flush << endl;
        }

        dash::finalize();
        exit(0);
    }
}

void initgrid( MatrixT& grid ) {

    if ( 0 != dash::myid() ) return;

    /* do initialization only from unit 0 with global accesses.
    This is not the fastest way and a parallel routine would not be too
    complicated but it should demonstrate the possibility to trade
    convenience against performance. */

    size_t d= grid.extent(2);
    size_t w= grid.extent(1);
    size_t h= grid.extent(0);

    size_t sum = 0;

    cout << "assign values to finest level " << w << " x " << h << " x " << d << endl;

    for( size_t j = 1; j < d-1; ++j ) {
    for( size_t i = 1; i < h-1; ++i ) {
    for( size_t k = 1; k < w-1; ++k ) {
        grid[i][k][j] = 7.0;
        sum++;
        cout << d*w*h << " / " << sum <<"\r" << std::flush;
    }
    }
    }

    // cold block inside of cube
    for( size_t j = 0; j < d; ++j ) {
    for( size_t i = h/5; i < h*3/4; ++i ) {
    for( size_t k = w*1/5; k < w*1/3; ++k ) {
        grid[i][k][j] = 3.0;
    }
    }
    }
    cout << "initgrid finished" << endl;
}

void setboundary( MatrixT& grid ) {

    if ( 0 != dash::myid() ) return;

    /* do initialization only from unit 0 with global accesses.
    This is not the fastest way and a parallel routine would not be too
    complicated but it should demonstrate the possibility to trade
    convenience against performance. */

    size_t d= grid.extent(2);
    size_t w= grid.extent(1);
    size_t h= grid.extent(0);

    /* set entire boundary to a low value */
    for ( size_t z= 1; z < d; ++z ) {
        for ( size_t y= 0; y < h; ++y ) {
            for ( size_t x= 0; x < w; ++x ) {
                grid[y  ][x  ][0  ]= 5.0;     // front plane
                grid[y  ][x  ][d-1]= 5.0;     // back plane
                grid[0  ][x  ][z  ]= 5.0;     // up plane
                grid[h-1][x  ][z  ]= 5.0;     // bottom plane
                grid[y  ][0  ][z  ]= 5.0;     // leftside plane
                grid[y  ][w-1][z  ]= 5.0;     // rightside plane
            }
        }
    }

    /* set parts of boundary to a high value */
    for ( size_t z= 0; z < d; ++z ) {
        for ( size_t y= 0; y < h*1/5; ++y ) {
            for ( size_t x= 0; x < w/3; ++x ) {
                grid[y][x][0]= 10.0;
                grid[0][x][z]= 10.0;
                grid[y][0][z]= 10.0;
            }
        }
    }

    /* set other parts of boundary to a very low value */
    /* upper part */
    for ( size_t z= 0; z < d*1/5; ++z ) {
        for ( size_t y= h*1/5; y < h*2/5; ++y ) {
            for ( size_t x= w*4/5; x < w-1; ++x ) {
                grid[y][x  ][0  ]= 0.0;     // front
                grid[y][w-1][z]= 0.0;   // right side
            }
        }
    }
    /* lower part */
    for ( size_t z= 0; z < d*1/5; ++z ) {
        for ( size_t y= h*3/5; y < h*4/5; ++y ) {
            for ( size_t x= w*4/5; x < w-1; ++x ) {
                grid[y][x  ][0  ]= 0.0;     // front
                grid[y][w-1][z]= 0.0;   // right side
            }
        }
    }
}

void scaledownboundary( const MatrixT& fine, MatrixT& coarse ) {

    assert( (coarse.extent(1)-2)*2 +2 == fine.extent(2) );
    assert( (coarse.extent(1)-2)*2 +2 == fine.extent(1) );
    assert( (coarse.extent(0)-2)*2 +2 == fine.extent(0) );

    size_t dc= coarse.local.extent(2);
    size_t wc= coarse.local.extent(1);
    size_t hc= coarse.local.extent(0);
    size_t df= fine.local.extent(2);
    size_t wf= fine.local.extent(1);
    size_t hf= fine.local.extent(0);
    std::array< long int, 3 > cornerc= coarse.pattern().global( {0,0,0} );

    size_t startx= ( 0 == cornerc[1] % 2 ) ? 1 : 0;
    size_t starty= ( 0 == cornerc[0] % 2 ) ? 1 : 0;
    size_t startz= ( 0 == cornerc[2] % 2 ) ? 1 : 0;

    if ( cornerc[2] == 0 )
    {
        // front plane
        for ( size_t y = 0; y < hc-starty; ++y )
        {
            for ( size_t x = 0; x < wc-startx; ++x )
            {
                coarse.local[starty+y][startx+x][0] = 0.25 *
                   (fine.local[starty+y*2  ][startx+x*2  ][0] +
                    fine.local[starty+y*2+1][startx+x*2  ][0] +
                    fine.local[starty+y*2  ][startx+x*2+1][0] +
                    fine.local[starty+y*2+1][startx+x*2+1][0]);
            }
        }
    }
    else // cornerc[2] != 0
    {
        // back plane -- only true for case with 8 units
        for ( size_t y = 0; y < hc-starty; ++y )
        {
            for ( size_t x = 0; x < wc-startx; ++x )
            {
                coarse.local[starty+y][startx+x][dc-1] = 0.25 *
                   (fine.local[starty+y*2  ][startx+x*2  ][df-1] +
                    fine.local[starty+y*2+1][startx+x*2  ][df-1] +
                    fine.local[starty+y*2  ][startx+x*2+1][df-1] +
                    fine.local[starty+y*2+1][startx+x*2+1][df-1]);
            }
        }
    }

    if ( cornerc[1] == 0 )
    {
        // left side
        for ( size_t z = 0; z < dc-startz; ++z )
        {
            for ( size_t y = 0; y < hc-starty; ++y )
            {
                coarse.local[starty+y][0][startz+z] = 0.25 *
                   (fine.local[starty+y*2  ][0][startz+z*2  ] +
                    fine.local[starty+y*2+1][0][startz+z*2  ] +
                    fine.local[starty+y*2  ][0][startz+z*2+1] +
                    fine.local[starty+y*2+1][0][startz+z*2+1]);
            }
        }
    }
    else
    {
        // right side -- only true for case with 8 units
        for ( size_t z = 0; z < dc-startz; ++z )
        {
            for ( size_t y = 0; y < hc-starty; ++y )
            {
                coarse.local[starty+y][wc-1][startz+z] = 0.25 *
                   (fine.local[starty+y*2  ][wf-1][startz+z*2  ] +
                    fine.local[starty+y*2+1][wf-1][startz+z*2  ] +
                    fine.local[starty+y*2  ][wf-1][startz+z*2+1] +
                    fine.local[starty+y*2+1][wf-1][startz+z*2+1]);
            }
        }
    }

    if ( cornerc[0] == 0 )
    {
        // top
        for ( size_t z = 0; z < dc-startz; ++z )
        {
            for ( size_t x = 0; x < wc-startx; ++x )
            {
                coarse.local[0][startx+x][startz+z] = 0.25 *
                   (fine.local[0][startx+x*2  ][startz+z*2  ] +
                    fine.local[0][startx+x*2+1][startz+z*2  ] +
                    fine.local[0][startx+x*2  ][startz+z*2+1] +
                    fine.local[0][startx+x*2+1][startz+z*2+1]);
            }
        }
    }
    else
    {
        // bottom -- only true in case with 8 units
        for ( size_t z = 0; z < dc-startz; ++z )
        {
            for ( size_t x = 0; x < wc-startx; ++x )
            {
                coarse.local[hc-1][startx+x][startz+z] = 0.25 *
                   (fine.local[hf-1][startx+x*2  ][startz+z*2  ] +
                    fine.local[hf-1][startx+x*2+1][startz+z*2  ] +
                    fine.local[hf-1][startx+x*2  ][startz+z*2+1] +
                    fine.local[hf-1][startx+x*2+1][startz+z*2+1]);
            }
        }
    }
}

void scaledown( Level& fine, Level& coarse ) {

    assert( (coarse.grid.extent(1)-2)*2 +2 == fine.grid.extent(1) );
    assert( (coarse.grid.extent(0)-2)*2 +2 == fine.grid.extent(0) );
    assert( (coarse.grid.extent(2)-2)*2 +2 == fine.grid.extent(2) );

    fine.halo.updateHalosAsync();

    const std::array< long unsigned int, 3 > extentc= coarse.grid.pattern().local_extents();
    const std::array< long signed int, 3 >   cornerc= coarse.grid.pattern().global( {0,0,0} );
    const std::array< long unsigned int, 3 > extentf= fine.grid.pattern().local_extents();

    size_t startx= ( 0 == cornerc[1] % 2 ) ? 1 : 0;
    size_t starty= ( 0 == cornerc[0] % 2 ) ? 1 : 0;
    size_t startz= ( 0 == cornerc[2] % 2 ) ? 1 : 0;


    for ( size_t z= 0; z < extentc[2]-startz; ++z ) {
        for ( size_t y= 0; y < extentc[0]-starty; ++y ) {
            for ( size_t x= 0; x < extentc[1]-startx; ++x ) {

                coarse.grid.local[starty+y][startx+x][startz+z]= 0.125 * (
                    fine.grid.local[starty+2*y  ][startx+2*x  ][startz+2*z  ] +
                    fine.grid.local[starty+2*y  ][startx+2*x+1][startz+2*z  ] +
                    fine.grid.local[starty+2*y+1][startx+2*x  ][startz+2*z  ] +
                    fine.grid.local[starty+2*y  ][startx+2*x  ][startz+2*z+1] +
                    fine.grid.local[starty+2*y+1][startx+2*x  ][startz+2*z+1] +
                    fine.grid.local[starty+2*y  ][startx+2*x+1][startz+2*z+1] +
                    fine.grid.local[starty+2*y+1][startx+2*x+1][startz+2*z  ] +
                    fine.grid.local[starty+2*y+1][startx+2*x+1][startz+2*z+1] );
            }
        }
    }

    fine.halo.waitHalosAsync();


    if ( 0 == cornerc[0] ) {

        // cout << "first row globally is boundary condition, nothing to do here" << endl;

    } else if ( 0 == cornerc[0] % 2 ) {

        cout << dash::myid() << ": first row locally is even, thus need to do scale down with halo row" << endl;

    } // else nothing to do here


    if ( coarse.grid.extent(0) == cornerc[0] + extentc[0] ) {

        // cout << "last row globally is boundary condition, nothing to do here" << endl;

    } else if ( 0 == ( cornerc[0] + extentc[0] ) % 2 ) {

        cout << dash::myid() << ": last row locally is odd (following first row in next block is even), thus need to do scale down with halo row" << endl;
    } // else nothing to do

    if ( 0 == cornerc[1] ) {

        // cout << "first column globally is boundary condition, nothing to do here" << endl;

    } else if ( 0 == cornerc[1] % 2 ) {

        cout << dash::myid() << ": first column locally is even, thus need to do scale down with halo column" << endl;

    } // else nothing to do here


    if ( coarse.grid.extent(1) == cornerc[1] + extentc[1] ) {

        // cout << "last column globally is boundary condition, nothing to do here" << endl;

    } else if ( 1 == ( cornerc[1] + extentc[1] ) % 2 ) {

        cout << dash::myid() << ": last column locally is odd (following first column in next block is even), thus need to do scale down with halo column" << endl;
    } // else nothing to do



    /* now do border elements if necessary */

    /* change to inner/outer scheme
         - nonblocking halo exchange if needed
         - check local indices to find start at even positions per dimension
         - check if first row/column needs to use halo
         - check if last row/column needs to use halo

    */


    dash::barrier();
}


void scaleup( Level& coarse, Level& fine ) {

    assert( (coarse.grid.extent(1)-2)*2 +2 == fine.grid.extent(1) );
    assert( (coarse.grid.extent(0)-2)*2 +2 == fine.grid.extent(0) );
    assert( (coarse.grid.extent(2)-2)*2 +2 == fine.grid.extent(2) );

    size_t w= coarse.grid.local.extent(1);
    size_t h= coarse.grid.local.extent(0);
    size_t d= coarse.grid.local.extent(2);
    std::array< long int, 3 > cornerc= coarse.grid.pattern().global( {0,0,0} );
    size_t startx= ( 0 == cornerc[1] ) ? 1 : 0;
    size_t starty= ( 0 == cornerc[0] ) ? 1 : 0;
    size_t startz= ( 0 == cornerc[2] ) ? 1 : 0;

    for ( size_t z= 0; z < d-1; ++z ) {
        for ( size_t y= 0; y < h-1 ; ++y ) {
            for ( size_t x= 0; x < w-1 ; ++x ) {

                double t= coarse.grid.local[starty+y][startx+x][startz+z];
                fine.grid.local[starty+2*y  ][startx+2*x  ][startz+2*z  ]= t;
                fine.grid.local[starty+2*y  ][startx+2*x  ][startz+2*z+1]= t;
                fine.grid.local[starty+2*y  ][startx+2*x+1][startz+2*z  ]= t;
                fine.grid.local[starty+2*y+1][startx+2*x  ][startz+2*z  ]= t;
                fine.grid.local[starty+2*y+1][startx+2*x  ][startz+2*z+1]= t;
                fine.grid.local[starty+2*y+1][startx+2*x+1][startz+2*z  ]= t;
                fine.grid.local[starty+2*y  ][startx+2*x+1][startz+2*z+1]= t;
                fine.grid.local[starty+2*y+1][startx+2*x+1][startz+2*z+1]= t;
            }
        }
    }

    dash::barrier();
}

double smoothen( Level& level ) {

    double res= 0.0;

    /// relaxation coeff.
    const double c= 1.0;

    /// begin pointer of local block, needed because halo border iterator is read-only
    auto gridlocalbegin= level.grid.lbegin();

    // async halo update
    level.halo.updateHalosAsync();

    // update inner

    size_t lw= level.grid.local.extent(1);
    size_t lh= level.grid.local.extent(0);
    size_t ld= level.grid.local.extent(2);

    /* the start value for both, the y loop and the x loop is 1 because either there is
    a border area next to the halo -- then the first column or row is covered below in
    the border update -- or there is an outside border -- then the first column or row
    contains the boundary values. */
    cout << "Smoothen start" << endl;
    for ( size_t z= 1; z < ld-1; z++ ) {
        for ( size_t y= 1; y < lh-1; y++ ) {

            // pointer movement manually
            /*const size_t d= 1;
            double* p_here= &gridlocalbegin[ d+ld*(y+lw*x) ];
            double* p_east= &gridlocalbegin[ d+ld*(y+lw*(x+1)) ];
            double* p_west= &gridlocalbegin[ d+ld*(y+lw*(x-1)) ];
            double* p_north= &gridlocalbegin[ d+ld*(y-1+lw*x) ];
            double* p_south= &gridlocalbegin[ d+ld*(y+1+lw*x) ];
            double* p_far= &gridlocalbegin[ d+1+ld*(y+lw*x) ];
            double* p_close= &gridlocalbegin[ d-1+ld*(y+1+lw*x) ];*/

            for ( size_t x= 1; x < lw-1; x++ ) {

                double dtheta =
                  ( level.grid.local[y  ][x+1][z  ] + // east
                    level.grid.local[y+1][x  ][z  ] +   // south
                    level.grid.local[y  ][x  ][z+1] +   // far
                    level.grid.local[y-1][x  ][z  ] +   // north
                    level.grid.local[y  ][x-1][z  ] +   // west
                    level.grid.local[y  ][x  ][z-1] )   // close
                    / 6.0 - level.grid.local[y  ][x  ][z  ]; // centre

                level.grid.local[y][x][z] = level.grid.local[y][x][z]+( c * dtheta);
                res= std::max( res, std::fabs( dtheta ) );

                /*double dtheta = ( *p_east + *p_west + *p_north + *p_south + *p_far + *p_close ) / 6.0 - *p_here ;
                *p_here += c * dtheta;
                res= std::max( res, std::fabs( dtheta ) );
                p_here++;
                p_east++;
                p_west++;
                p_north++;
                p_south++;
                p_far++;
                p_close++;*/
            }
        }
    }
    cout << "smoothen finish" << endl;
    // wait for async halo update
    level.halo.waitHalosAsync();

    auto bend = level.halo.bend();
    // update border area
    for( auto it = level.halo.bbegin(); it != bend; ++it ) {

        double dtheta = ( it.valueAt(2) + it.valueAt(3) +
                          it.valueAt(0) + it.valueAt(1) +
                          it.valueAt(4) + it.valueAt(5) ) / 6.0 - *it;
        //gridlocalbegin[ it.lpos() ] += c * dtheta;
        *it += c*dtheta;
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

void v_cycle( vector<Level*>& levels, uint32_t inneriterations= 1 ) {

    writeToCsv( levels[0]->grid ); // 1st

    for ( auto i= 1; i < levels.size(); i++ ) {

        double res= smoothen( *levels[i-1] );

        writeToCsv( levels[i-1]->grid ); // mod 2th

        scaledown( *levels[i-1], *levels[i] );

        if ( 0 == dash::myid() ) {
            cout << "smoothen with residual " << res <<
                ", then scale down " << levels[i-1]->grid.extent(1) << "x" << levels[i-1]->grid.extent(1) <<
                " --> " << levels[i]->grid.extent(1)   << "x" << levels[i]->grid.extent(1) << endl;
        }

        writeToCsv( levels[i]->grid ); // mod 3rd

        //exit(1);
    }

    /* do a fixed number of smoothing steps until the residual on the coarsest grid is
    what you want on the finest level in the end. */
    for ( auto i= 0; i < inneriterations; i++ ) {

        double res= smoothen( *levels.back() );

        if ( 0 == dash::myid() ) {
            cout << "smoothen coarsest with residual " << res << endl;
        }
        writeToCsv( levels.back()->grid );
    }

    for ( auto i= levels.size()-1; i > 0; i-- ) {

        scaleup( *levels[i], *levels[i-1] );

        writeToCsv( levels[i-1]->grid );

        double res= smoothen( *levels[i-1] );

        if ( 0 == dash::myid() ) {
            cout << "scale up " << levels[i]->grid.extent(1)   << "x" << levels[i]->grid.extent(1) <<
                " --> " << levels[i-1]->grid.extent(1) << "x" << levels[i-1]->grid.extent(1) <<
                ", then smoothen with residual " << res << endl;
        }

        writeToCsv( levels[i-1]->grid );
    }
}

int main( int argc, char* argv[] )
{

    dash::init(&argc, &argv);

    TeamSpecT teamspec( dash::Team::All().size(), 1, 1);
    teamspec.balance_extents();

    constexpr uint32_t howmanylevels= 6;
    vector<Level*> levels;
    levels.reserve( howmanylevels );

    /* create all grid levels, starting with the finest and ending with 5 x 5 x X */
    for ( auto l= 0; l < howmanylevels-2; l++ ) {

        if ( 0 == dash::myid() ) {
            cout << "level " << l << " is " << (1<<(howmanylevels-l))+2 << "x" <<
            (1<<(howmanylevels-l))+2 << "x" << (1<<(howmanylevels-l))+2 << endl;
        }

        levels.push_back( new Level( (1<<(howmanylevels-l)), (1<<(howmanylevels-l)), (1<<(howmanylevels-l)),
            dash::Team::All(), teamspec ) );

        dash::barrier();
        cout << "unit " << dash::myid() << " : " <<
            levels.back()->grid.local.extent(1) << " x " <<
            levels.back()->grid.local.extent(0) << " x " <<
            levels.back()->grid.local.extent(2)<< endl;
        dash::barrier();


        if ( 0 == l ) {
            sanitycheck( levels[0]->grid );
            setboundary( levels[0]->grid );
            continue;
        }

        scaledownboundary( levels[l-1]->grid, levels[l]->grid );
    }

    dash::barrier();

    /* Fill finest level. Strictly, we don't need to set any initial values here
    but we do it for demonstration in the CSV files */
    initgrid( levels[0]->grid );

    dash::barrier();

    v_cycle( levels, 10 );

    double res= smoothen( *levels.front() );
    if ( 0 == dash::myid() ) {
        cout << "final residual " << res << endl;
    }

    writeToCsv( levels.front()->grid );

    dash::finalize();
}
