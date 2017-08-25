using dash::transform;
using dash::NArray;
using dash::Array;
using dash::max_element;
using dash::BLOCKED;
using dash::size;
using dash::barrier;
using dash::Shared;

template<typename T = MATRIX_T>
inline void Thresh(
  NArray< T   , 2 > const & rand_mat   ,
	NArray< bool, 2 >       & thresh_mask,
               uint const   nrows      ,
               uint const   ncols      ,
               uint const   percent    )
{
  // find max value in rand_mat
  auto max_glob = max_element( rand_mat.begin( ), rand_mat.end( ) );

  T max = *max_glob;

  // get number of units running
  size_t num_units = size( );

  // create global histo array and initialze with 0
  Array<uint> histo( (max + 1) * num_units, BLOCKED );

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
  barrier( );

  // add the values of the local histogram to the histogram of unit0
  if( 0 != myid ) {
    transform<uint>  ( histo.lbegin( ),
    histo.lend       ( )              ,
    histo.begin      ( )              , // points to global begin -> lbegin of unit0
    histo.begin      ( )              ,
    dash::plus<uint> ( )              );
  }

  // create new shared variable
  Shared<int> threshold;

  // wait for all units to finish adding (especially unit0 should wait)
  barrier( );

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
        if( histo.local[j] ) cout << std::setw(3) << j
        << " counted: " << histo.local[j] << endl;
      }
    #endif

    // count defines how many values are to be hold on given percentage
    uint count = ( static_cast<size_t>(nrows) * ncols * percent ) / 100;
    uint prefixsum = 0;
    int  i;

    // find threshold
    for( i = max; i >= 0 && prefixsum <= count; --i ) {
      prefixsum += histo.local[i];
    }

    threshold.set( ++i );
    #if 0
      cout << "original threshold: " << i << " - perc: " << percent << endl;
    #endif
  }

  // wait for unit0 to finish and flush
  threshold.barrier( );

  int threshLclCpy = threshold.get( );

  T const * src  = rand_mat.lbegin( );
  bool * i = thresh_mask.lbegin   ( );

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
  barrier( );
}