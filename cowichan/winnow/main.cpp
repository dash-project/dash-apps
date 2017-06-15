/* winnow: weighted point selection
 *
 * input:
 *   matrix: an integer matrix, whose values are used as masses
 *   mask: a boolean matrix, showing which points are eligible for
 *     consideration
 *   nrows, ncols: number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   points: a vector of (x, y) points
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <vector>

#include <libdash.h>
#define MATRIX_T uchar

using uint  = unsigned int ;
using uchar = unsigned char;
using Point = struct{ uint row, col;};
using namespace std;



struct InputPar { uint nrows, ncols; } in;
uint nelts;
static int myid;

#include "winnow_placeholder.h"


inline void ReadRowsNCols(){

  dash::Shared<InputPar> input_transfer;

  if(0 == myid)
  {
    cin >> in.nrows;
    cin >> in.ncols;

    input_transfer.set(in);
  }
  input_transfer.barrier();
  in = input_transfer.get();
}


template< typename T, typename X>
inline void ReadMatricesAndNelts( T& randMat, X& threshMask ){

  dash::Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    //read matrices
    MATRIX_T tmp;
      for ( auto i : randMat ){
        scanf( "%u", &tmp )  , i = tmp;
      }
      bool tmpB;
      for ( auto i : threshMask ){
        scanf( "%u", &tmpB ) , i = tmpB;
      }
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


int main(int argc, char** argv) {
  
  dash::init( &argc,&argv );
  
  myid = dash::myid( );
  ReadRowsNCols();
  
  dash::NArray< MATRIX_T, 2 > randMat    ( in.nrows, in.ncols );
  dash::NArray< bool,     2 > threshMask ( in.nrows, in.ncols );
  
  ReadMatricesAndNelts( randMat, threshMask );
  
  dash::Array< Point > points(nelts);
  
  // matrix = new vector<vector <int> >(nrows, vector<int>(ncols));
  // mask = new vector<vector<int> >(nrows, vector<int>(ncols));

  // read_matrix(nrows, ncols, matrix);
  // read_matrix(nrows, ncols, mask);

  // scanf("%d", &nelts);

  // vector<pair<int, int> > points(nelts);

  if(0 == myid){
    winnow( in.nrows, in.ncols, randMat, threshMask, nelts, points );

    printf("%d\n", nelts);

    for ( auto i : points) {
      printf("%d %d\n", static_cast<Point>(i).row, static_cast<Point>(i).col);
    }
    printf("\n");
  }
  dash::finalize( );
}
