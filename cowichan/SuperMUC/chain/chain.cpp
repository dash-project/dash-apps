#include <libdash.h>
#include <iostream>

using std::cout;
using std::endl;
using std::cin;

using dash::Shared;
using dash::NArray;
using dash::Array;

using uint     = unsigned int  ;
using uchar    = unsigned char ;
using MATRIX_T = uint          ;

static int myid;

#include "../randmat/randmat.h"
#include "../thresh/thresh.h"
// #include "../winnow/winnow_placeholder.h"
#include "../winnow/winnow.h"
#include "../outer/outer.h"
#include "../product/product.h"

#include <time.h>
#include <stdio.h>
#include <cstring>
using std::strcmp;

struct InputPar { uint nRowsCols, seed, thresh, winnow_nelts; } in;


/* One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadPars( )
{
  Shared<InputPar> input_transfer;

  if(0 == myid)
  {
    cin >> in.nRowsCols    ;
    cin >> in.seed         ;
    cin >> in.thresh       ;
    cin >> in.winnow_nelts ;

    input_transfer.set(in) ;
  }
  input_transfer.barrier() ;
  in = input_transfer.get();
}

/* both parameters musst have same size!
 * copy from a dash::Array to std::vector
 */
inline void CopyFromDashToStd(
  Array  <double> const & dashVector,
  vector <double>       & loclVector)
{
  if(0==myid) dash::copy(dashVector.begin(), dashVector.end(), loclVector.data());
  
  BroadcastOuterVecToUnits( loclVector ); //defined in product.h
}

/* both parameters musst have same size!
 * copy from a dash::Array to std::vector
 */
template<typename T, typename X, typename Y>
inline void CopyFromDashToStd(
  Array<T, X, Y>  const & dashVector,
       vector<T>        & loclVector)
{
  if(0==myid) dash::copy(dashVector.begin(), dashVector.end(), loclVector.data());

  dart_ret_t ret = dart_bcast(
                    static_cast<void*>( loclVector.data() ),  // buf
                    loclVector.size() * sizeof(T)          ,  // nelem
                    DART_TYPE_BYTE                         ,  // dtype
                    dash::team_unit_t(0)                   ,  // root/source
                    dash::Team::All( ).dart_id( )          ); // team

  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl;
}


// @echo $(nRowsCols) $(seed) $(thresh) $(winnow_nelts) | $(TBB_ROOT)/chain/expertpar/main
int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );

  struct timespec start, stop;
  double accum;
  int is_bench = 0;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    }
  }

  myid = dash::myid( );
  ReadPars( );

  // initialize variables
  NArray<MATRIX_T, 2> rand_mat      ( in.nRowsCols   , in.nRowsCols    );
  NArray<bool    , 2> thresh_mask   ( in.nRowsCols   , in.nRowsCols    );
  NArray<double  , 2> outer_mat     ( in.winnow_nelts, in.winnow_nelts );
  Array <double     > outer_vec     ( in.winnow_nelts                  );
  vector<double     > prod_vec      ( in.winnow_nelts                  );
  // Array  < double             > result        ( in.winnow_nelts                  );
  res_array_t winnow_dashRes; // defined in winnow.h

  // after the run of outer "outer_vec" will be recycled/reused for the final output
  auto & result = outer_vec;
  
  vector<value> winnow_vecRes (in.winnow_nelts); //value defined in winnow.h

  
  if( clock_gettime( CLOCK_MONOTONIC, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  

  Randmat( rand_mat, in.seed );
  Thresh( rand_mat, thresh_mask, in.nRowsCols, in.nRowsCols, in.thresh );

  Winnow( in.nRowsCols, in.nRowsCols, rand_mat, thresh_mask, in.winnow_nelts, winnow_dashRes );
  
  CopyFromDashToStd( winnow_dashRes, winnow_vecRes );

  
  Outer( winnow_vecRes, outer_mat, outer_vec, in.winnow_nelts );
  
  CopyFromDashToStd( outer_vec, prod_vec );
  

  Product( prod_vec, outer_mat, result, in.winnow_nelts );
  if( clock_gettime( CLOCK_MONOTONIC, &stop) == -1 ) {
    perror( "clock gettime error 2" );
    exit( EXIT_FAILURE );
  }

  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;

  if( 0 == myid ){
    // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
    printf( "DASH  ,Chain  ,%5u,%5u,%3u,%5u,%2lu,%.9lf,isBench:%d\n", in.nRowsCols, in.nRowsCols, in.thresh, in.winnow_nelts, dash::Team::All().size(), accum, is_bench );
  }

  if( !is_bench ){ PrintOutput( result, in.winnow_nelts ); }
  dash::finalize( );
}