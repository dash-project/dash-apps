using dash::NArray;
using dash::Array;
using dash::max_element;
using dash::BLOCKED;
using dash::size;
using dash::barrier;
using dash::Shared;

template<typename T = MATRIX_T,typename hisT = size_t>
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
  Array<hisT> histo( (max + 1) * num_units, BLOCKED );

  // initialize the histogram
  for( hisT * i = histo.lbegin(); i < histo.lend(); ++i) {
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
    dash::transform(
      histo.lbegin     ( ) ,
      histo.lend       ( ) ,
      histo.begin      ( ) , // points to global begin -> lbegin of unit0
      histo.begin      ( ) ,
      dash::plus<hisT> ( ) );
  }

  // create new shared variable
  Shared<T> threshold;

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

    // count defines how many values are to be hold on given percentage
    hisT count = ( static_cast<size_t>(nrows) * ncols * percent ) / 100;
    hisT prefixsum = 0;
    int  i;

    // find threshold
    for( i = max; i >= 0 && prefixsum <= count; --i ) {
      prefixsum += histo.local[i];
    }

    threshold.set( ++i );
  }
  // wait for unit0 to finish and flush
  threshold.barrier( );

  T threshLclCpy = threshold.get( );

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

  // wait for all units finish calculating local boolean mask
  barrier( );
}