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

struct InputPar { uint nRowsCols, seed, thresh, winnow_nelts; } in;


/* 
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadPars(){

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


// @echo $(nRowsCols) $(seed) $(thresh) $(winnow_nelts) | $(TBB_ROOT)/chain/expertpar/main
int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );

  myid = dash::myid( );
  ReadPars( );

  // initialize variables
  NArray < MATRIX_T, 2        > rand_mat      ( in.nRowsCols, in.nRowsCols       );
  NArray < bool    , 2        > thresh_mask   ( in.nRowsCols, in.nRowsCols       );
  vector < pair<POI_T, POI_T> > winnow_points ( in.winnow_nelts                  );
  NArray < double  , 2        > outer_mat     ( in.winnow_nelts, in.winnow_nelts );
  Array  < double             > outer_vec     ( in.winnow_nelts                  );

  // execute functions
  Randmat( rand_mat, in.nRowsCols, in.nRowsCols, in.seed );
  Thresh ( rand_mat, thresh_mask, in.nRowsCols, in.nRowsCols, in.thresh );
  if( 0 == myid ){
    winnow( in.nRowsCols, in.nRowsCols, rand_mat, thresh_mask, in.winnow_nelts, winnow_points );
  }
  BroadcastPointsToUnits( winnow_points );
  Outer( winnow_points, outer_mat, outer_vec, in.winnow_nelts );

  dash::finalize( );
}