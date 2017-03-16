#include <libdash.h>
#include <iostream>
#include <cstdio>

using std::cout;
using std::endl;

typedef unsigned int  uint;
typedef unsigned char uchar;

#define MTRX_TYPE uchar


template< typename T = MTRX_TYPE >
inline void thresh(uint nrows, uint ncols, uint percent, int myid);


template< typename T >
void print2d( T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>( mat(i,j) )<< " ";
    }
    cout << endl;
  }
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc, &argv );
  int myid = static_cast<int>( dash::myid( ) );
  
  if( argc != 4 ){
    if( 0 == myid ){ cout << "3 Parameters expected!"         << endl
                    << "Usage: cowichan_thresh nRows nCols percentage" << endl
                    << "Then enter the matrix."               << endl;
    }
    dash::finalize( );
    return 0;
  }
  
  uint nrows   = static_cast<uint>( atoi( argv[1] ) );
  uint ncols   = static_cast<uint>( atoi( argv[2] ) );
  uint percent = static_cast<uint>( atoi( argv[3] ) );
  
  thresh( nrows, ncols, percent, myid );

  dash::finalize( );
}


template<typename T = MTRX_TYPE>
inline void thresh(uint nrows, uint ncols, uint percent, int myid){

  dash::NArray<T   , 2> matSrc( nrows, ncols );
  dash::NArray<bool, 2> mask  ( nrows, ncols );
  
  // read input matrix from stdin
  if( 0 == myid ){
    T tmp;
    for ( auto i : matSrc ){
      scanf( "%d", &tmp );
      i = tmp;
    }
  }
  
  // wait for initialization of the matrix before calculating maximum
  dash::barrier( );
  
  // find max value in matSrc
  auto maxGlobIt = dash::max_element( matSrc.begin( ), matSrc.end( ) );
  T max          = *maxGlobIt;
  T maxPO        =  max + 1  ;   //max plus one
  
  // get number of units running
  size_t num_units = dash::size( );
  
  // create global histo array and initialze with 0
  dash::Array<uint> histo( maxPO * num_units, dash::BLOCKED );
  for( uint * i = histo.lbegin(); i < histo.lend(); ++i) {
    *i = 0;
  }
  
  for( T * i = matSrc.lbegin( ); i < matSrc.lend( ); ++i ) {
    ++histo.local[*i];
  }
  
  /* barrier is necessary because if unit 0 is still calculating
     while another unit starts with dash::transform there could be a race condition */
  dash::barrier( );
  
  // add the values of the local histogram to the histogram of unit 0
  if( 0 != myid ) {
    dash::transform<uint>( histo.lbegin     ( ) ,
                           histo.lend       ( ) ,
                           histo.begin      ( ) ,
                           histo.begin      ( ) ,
                           dash::plus<uint> ( ));
  }
  
  // wait for all units to finish adding
  dash::barrier( );
  
  // create new shared variable
  dash::Shared<int> threshold;
  
  if( 0 == myid ) {
    
    #if 0
      for( uint j = 0; j < histo.lsize( ); ++j ) {
        if( histo.local[j] ) cout << std::setw(3)   << j 
                  << " counted: " << histo.local[j] << endl;
      }
    #endif

    uint count     = ( nrows * ncols * percent ) / 100;
    uint prefixsum = 0;
    int  i;

    // find threshold
    for( i = max; i >= 0 && prefixsum <= count; --i ) {
      prefixsum += histo.local[i];
    }
    
    threshold.set( ++i );
    #if 0
      cout << "original threshold: " << i << endl;
    #endif
  }
  
  // wait for Unit 0 to finish and flush
  threshold.barrier( );

  {
  int threshLclCpy = threshold.get( );
  
    T * src  = matSrc.lbegin( );
    bool * i = mask.lbegin(   );
    
        *i =   *src >= threshLclCpy;      
    while ( i < mask.lend( ) ) {
      *++i = *++src >= threshLclCpy;
    }
  #if 0
    cout << myid << " got threshold: " << threshLclCpy << endl;
  #endif
  }
  
  // wait for all units finish calculating local boolean mask
  dash::barrier( );
  if( 0 == myid ) {
    cout << "----------------------------" << endl;
    print2d(mask);
  }
}





















