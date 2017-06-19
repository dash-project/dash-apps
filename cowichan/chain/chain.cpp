#include <libdash.h>
#include <iostream>

using std::cout;
using std::endl;
using std::cin;

using dash::Shared;
using dash::NArray;
using dash::Array;

using uint  = unsigned int ;
using uchar = unsigned char;
using MATRIX_T = uchar;
using POI_T = int;  //this type musst be signed!

static int myid;

#include "../randmat/randmat.h"
#include "../thresh/thresh.h"
#include "../winnow/winnow_placeholder.h"
#include "../outer/outer.h"
#include "../product/product.h"

struct InputPar { uint nRowsCols, seed, thresh, winnow_nelts; } in;


/* 
 * One unit has the job to read in the parameters.
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
  BroadcastOuterVecToUnits( loclVector );
}


// @echo $(nRowsCols) $(seed) $(thresh) $(winnow_nelts) | $(TBB_ROOT)/chain/expertpar/main
int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );

  myid = dash::myid( );
  ReadPars( );

  // initialize variables
  NArray < MATRIX_T, 2        > rand_mat      ( in.nRowsCols   , in.nRowsCols    );
  NArray < bool    , 2        > thresh_mask   ( in.nRowsCols   , in.nRowsCols    );
  vector < pair<POI_T, POI_T> > winnow_points ( in.winnow_nelts                  );
  NArray < double  , 2        > outer_mat     ( in.winnow_nelts, in.winnow_nelts );
  Array  < double             > outer_vec     ( in.winnow_nelts                  );
  vector < double             > prod_vec      ( in.winnow_nelts                  );
  // Array  < double             > result        ( in.winnow_nelts                  );
  
  // after the run of outer "outer_vec" will be recycled/reused for the final output
  auto & result = outer_vec;
  
  
  // execute functions
  Randmat( rand_mat, in.seed );
  
  Thresh ( rand_mat, thresh_mask, in.nRowsCols, in.nRowsCols, in.thresh );
  
  if(0 == myid){ winnow( in.nRowsCols, in.nRowsCols, rand_mat, thresh_mask, in.winnow_nelts, winnow_points );}
  BroadcastPointsToUnits( winnow_points );
  
  Outer( winnow_points, outer_mat, outer_vec, in.winnow_nelts );
  CopyFromDashToStd( outer_vec, prod_vec );
  
  Product( prod_vec, outer_mat, result, in.winnow_nelts );
  
  PrintOutput( result, in.winnow_nelts );

  dash::finalize( );
}