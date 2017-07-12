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
using MATRIX_T = uchar         ;

static int myid;

#include "../randmat/randmat.h"
#include "../thresh/thresh.h"
// #include "../winnow/winnow_placeholder.h"
#include "../winnow/winnow.h"
#include "../outer/outer.h"
#include "../product/product.h"

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
  if(0 == myid)
  {
    double * vec = loclVector.data( );
    for( double const i : dashVector )
    {
      *(vec++) = i;
    }
  }
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
  if(0 == myid)
  {
    T * vec = loclVector.data( );
    for( T const i : dashVector )
    {
      *(vec++) = i;
    }
  }
  
  // using dash::Team        ;
  // using dash::team_unit_t ;
  
  // Team & team = dash::Team::All( );
                    
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

  myid = dash::myid( );
  ReadPars( );

  // initialize variables
  NArray<MATRIX_T, 2> rand_mat      ( in.nRowsCols   , in.nRowsCols    );
  NArray<bool    , 2> thresh_mask   ( in.nRowsCols   , in.nRowsCols    );
  vector<value      > winnow_vecRes ( in.winnow_nelts                  ); //value defined in winnow.h
  NArray<double  , 2> outer_mat     ( in.winnow_nelts, in.winnow_nelts );
  Array <double     > outer_vec     ( in.winnow_nelts                  );
  vector<double     > prod_vec      ( in.winnow_nelts                  );
  // Array  < double             > result        ( in.winnow_nelts                  );
  res_array_t winnow_dashRes; // defined in winnow.h
  
  // after the run of outer "outer_vec" will be recycled/reused for the final output
  auto & result = outer_vec;
  
  
  // execute functions
  Randmat( rand_mat, in.seed );
  
  Thresh( rand_mat, thresh_mask, in.nRowsCols, in.nRowsCols, in.thresh );
  
  // if(0 == myid){ winnow( in.nRowsCols, in.nRowsCols, rand_mat, thresh_mask, in.winnow_nelts, winnow_points );}
  // BroadcastPointsToUnits( winnow_points );
  Winnow( in.nRowsCols, in.nRowsCols, rand_mat, thresh_mask, in.winnow_nelts, winnow_dashRes );
  CopyFromDashToStd( winnow_dashRes, winnow_vecRes );
  
  Outer( winnow_vecRes, outer_mat, outer_vec, in.winnow_nelts );
  CopyFromDashToStd( outer_vec, prod_vec );
  
  Product( prod_vec, outer_mat, result, in.winnow_nelts );
  
  PrintOutput( result, in.winnow_nelts );

  dash::finalize( );
}