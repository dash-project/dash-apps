#include <libdash.h>
#include <iostream>
#include <cstdio>

using std::cout;
using std::endl;
using std::cin;

using uint  = unsigned int ;
using uchar = unsigned char;

struct InputPar { uint nrows, ncols;} in;
uint percent;

static int myid;

#define MATRIX_T int

/*
 * This function prints the content of a 2D matrix to std::out.
 * Datatypes are casted to <const uint> for readable output
 * (otherwise uchars would be printed as chars and not as numerics)
 */
template<typename T>
inline void Print2D(const T& mat ) {
  if(0==myid){

    //cout << in.nrows << " " << in.ncols << "\n";

    for( int i = 0; i < mat.extent(0); i++ ) {
      for( int j = 0; j < mat.extent(1); j++ ) {
	cout << static_cast<const uint>( mat(i,j) )<< " ";
      }
      cout << "\n";
    }
  }
}

// unit0 reads the matrix with random values from std::in
template<typename T>
inline void ReadRandMat(dash::NArray<T,2>& rand_mat){
  if(0 == myid){
    T tmp;
    for ( auto i : rand_mat ){
      scanf( "%d", &tmp );
      i = tmp;
    }
  }

  // wait for initialization of the matrix before continuing
  dash::barrier( );
}

/*
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
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

inline void ReadPercentage(){

  dash::Shared<uint> percent_transfer;

  if(0 == myid)
  {
    cin >> percent;

    percent_transfer.set(percent);
  }
  percent_transfer.barrier();
  percent = percent_transfer.get();
}


template<typename T = MATRIX_T>
inline void Thresh(
  const dash::NArray<T   ,2>& rand_mat,
	dash::NArray<bool,2>& thresh_mask,
  const uint& nrows,
  const uint& ncols,
  const uint& percent
){
  // find max value in rand_mat
  auto max_glob = dash::max_element( rand_mat.begin( ), rand_mat.end( ) );

  T max = *max_glob;

  // get number of units running
  size_t num_units = dash::size( );

  // create global histo array and initialze with 0
  dash::Array<uint> histo( (max + 1) * num_units, dash::BLOCKED );

  // initialize the histogram
  for( uint * i = histo.lbegin(); i < histo.lend(); ++i) {
    *i = 0;
  }

  // every unit generates a histogram for the local values
  for( T const * i = rand_mat.lbegin( ); i < rand_mat.lend( ); ++i ) {
    ++histo.local[*i];
  }

  /* barrier is necessary because if unit 0 is still calculating
   * while another unit starts with dash::transform there could be a race condition
   */
  dash::barrier( );

  // add the values of the local histogram to the histogram of unit0
  if( 0 != myid ) {
    dash::transform<uint>( histo.lbegin     ( ) ,
			   histo.lend       ( ) ,
			   histo.begin      ( ) , // points to global begin -> lbegin of unit0
			   histo.begin      ( ) ,
			   dash::plus<uint> ( ));
  }

  // create new shared variable
  dash::Shared<int> threshold;

  // wait for all units to finish adding (especially unit0 should wait)
  dash::barrier( );

/*
 * In the following scope unit0 calculates the threshold for the
 * matrix with random values. A given percentage defines how much values
 * are to be hold in the result. Lower values are dropped first and so
 * unit0 calculates which "low" values are dropped and defines therefore
 * a threshold.
 */
  if( 0 == myid ) {

    // if compiled; unit0 prints the global histogram (which resides at this point only on unit0)
    #if 0
      for( uint j = 0; j < histo.lsize( ); ++j ) {
	if( histo.local[j] ) cout << std::setw(3)   << j
		  << " counted: " << histo.local[j] << endl;
      }
    #endif

    // count defines how many values are to be hold on given percentage
    uint count     = ( nrows * ncols * percent ) / 100;
    uint prefixsum = 0;
    int  i;

    // find threshold
    for( i = max; i >= 0 && prefixsum <= count; --i ) {
      prefixsum += histo.local[i];
    }

    threshold.set( ++i );
    #if 0
      cout << "original threshold: " << i << "perc: " << percent << endl;
    #endif
  }

  // wait for unit0 to finish and flush
  threshold.barrier( );

  int threshLclCpy = threshold.get( );

  T const * src  = rand_mat.lbegin( );
  bool * i = thresh_mask.lbegin(   );

  /* //debug ausgabe von rand_mat
  int co = 0;
  while( src < rand_mat.lend()){
    cout << *src++ << " ";
    if (++co == 9){
     cout << endl;
     co = 0;
    }
  }*/
  
  *i =  ( *(src) >= threshLclCpy);
  while ( i < thresh_mask.lend( ) - 1 ) {
    *(++i) = (*(++src) >= threshLclCpy);
  }

  #if 0
    cout << myid << " got threshold: " << threshLclCpy << endl;
  #endif

  // wait for all units finish calculating local boolean mask
  dash::barrier( );
}


int main( int argc, char* argv[] )
{
  dash::init( &argc, &argv );

  myid = dash::myid( );
  ReadRowsNCols( );

  dash::NArray<MATRIX_T,2> rand_mat   ( in.nrows, in.ncols );
  dash::NArray<bool    ,2> thresh_mask( in.nrows, in.ncols );

  ReadRandMat(rand_mat);
  ReadPercentage();

  Thresh(rand_mat, thresh_mask, in.nrows, in.ncols, percent);

  Print2D( thresh_mask );

  dash::finalize( );
}
