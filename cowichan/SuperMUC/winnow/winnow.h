#include "../Terminal_Color.h"

#include <cstdlib>  // for std::malloc
#include <cstring>  // for std::memcpy and std::memset

// #define DEBUG
// #define DEBUG_DETAILED  // requires DEBUG


#ifdef DEBUG
  #include <chrono>
  #include <thread>
#endif


#define MAX_KEY               999
#define MIN_NUM_ELEM_PER_UNIT 10


using std::cout         ;
using std::cin          ;
using std::endl         ;
using std::vector       ;
using std::pair         ;

using dash::Team        ;
using dash::Array       ;
using dash::NArray      ;
using dash::Shared      ;
using dash::fill        ;
using dash::team_unit_t ;

using uint        = unsigned  int ;
using uchar       = unsigned char ;
using MATRIX_T    =          uint ;
using HISTOG_T    =          long ;

using value       = struct{ uint row, col; }; //hast to be signed! - results are stored in here

using Point       = struct{ MATRIX_T value; uint row, col;   };  // results are created using these
using unitRange   = struct{ MATRIX_T begin, end; HISTOG_T count; };  // unitValueRange

using pattern_t   = dash::CSRPattern< 1, dash::ROW_MAJOR, long >;
using extent_t    = pattern_t::size_type;
using res_array_t = dash::Array<value, long, pattern_t>;


bool operator<(const Point& lhs, const Point& rhs)
{
  return lhs.value < rhs.value;
}


std::ostream& operator<<(std::ostream& os, const Point& p)
{
  #ifdef DEBUG
    // return os << "(" << fmt( p.value, FCYN ) << "," << fmt( p.row, FGREEN ) << "," << fmt( p.col, FGREEN ) << ")-";
    return os << "(" << fmt( p.value, FCYN ) << "," << fmt( p.row, FGREEN ) << "," << fmt( p.col, FGREEN ) << ")";
  #else
    return os << p.row << " " << p.col << "\n";
  #endif
}

std::ostream& operator<<(std::ostream& os, const value& p)
{
    return os << p.row << " " << p.col << "\n";
}


// #ifdef DEBUG
  using std::this_thread::sleep_for;
  inline void __sleep( uint const baseDur = 0, uint const mult = 10 )
  {
     uint SLEEP_TIME__ = (myid + 1) * mult + baseDur;
      sleep_for(std::chrono::milliseconds(SLEEP_TIME__));
  }
// #endif


template<typename T = MATRIX_T, typename hisT = HISTOG_T >
inline void Winnow(
              uint const   nrows      ,
              uint const   ncols      ,
          NArray<   T,2> & randMat    ,
          NArray<bool,2> & threshMask ,
              uint const   nelts      ,
       res_array_t       & result     )
{
  Team & team   = dash::Team::All ( );
  size_t nUnits = team.size       ( );


  /* create global histo array for sorting
   * size += 1 for direct Index access -> histo[2]++ counts for value 2
   * size += 1 for additional value at the end of the
   * histogram -> used for counter how many values were found
   */
  Array<hisT> histo( (MAX_KEY + 2) * nUnits, dash::BLOCKED );
  // fill( histo.begin( ), histo.end( ), 0 );
  // for(uint * it = histo.lbegin( ); it < histo.lend( ); ++it){*it = 0;}
  std::memset( histo.lbegin(), 0, histo.lsize() * sizeof(hisT) );

  // local found points are gathered in this vector
  vector<Point> pointsLocal;


  // returns a object with the global row and column of the the local coordinates {0,0}
  auto globIndex = randMat.pattern( ).global( {0,0} );

  uint         gRow   = globIndex[0]           ;
  uint         gCol   = globIndex[1]           ;
  hisT       & found  = histo.local[MAX_KEY+1] ;
     T const * matrEl = randMat.lbegin( )      ;


 /* read in local part of mask - matrix combination
  * and while doing that:
  * generate histogram and count the values (-> found)
  */
  for ( bool const * maskEl = threshMask.lbegin( );  maskEl < threshMask.lend( );  ++maskEl , ++matrEl, ++gCol )
  {
    if( gCol == ncols ) gCol = 0, ++gRow; // end of row -> next row in matrix/mask
    if( *maskEl       )
    {
      pointsLocal.push_back( Point{ *matrEl, gRow, gCol } );
      ++histo.local[*matrEl];
    }
  }
  
  found = pointsLocal.size();

  // is here a barrier needed? -> to wait for completion on unit0?! -> YES YES YES!
  dash::barrier();
  
  if( 0 != myid ) {
    dash::transform(
      histo.lbegin     ( ) ,
      histo.lend       ( ) ,
      histo.begin      ( ) ,
      histo.begin      ( ) ,
      dash::plus<hisT> ( ) );
  }

  // free randMat and threshMask
  randMat.deallocate();
  threshMask.deallocate();

  // in this array will be the distribution info for creating buckets
  unitRange * distr = static_cast<unitRange*>(  std::malloc( nUnits * sizeof(unitRange) )  );
  unitRange * const distr_end = distr + nUnits;

  
  // unit 0 have to wait for rest to finish adding their values
  dash::barrier();

  if( 0 == myid ) { // calculate bucket distribution
  
    /* calculate how much elements each unit should hold ideally
     * increment by one for safety garuantees in distribution
     * that is the assumption (ideal * nUnits > foundAllSize) == true
     */
     hisT ideal = std::max( static_cast<size_t>( MIN_NUM_ELEM_PER_UNIT ), (found / nUnits) + 1 );

    // begin - 1 for loop logic (starting with prefix ++)
    unitRange * uRPtr = distr;
    uRPtr->begin      = 0;

   /* the loop for calculation of the distribution got a bit more complex
    * because i wanted to iterate only once over "distr"
    * therefore no initialization beforehand is needed
    */
    if ( 1 == nUnits )
    {
      uRPtr->end   = MAX_KEY;
      uRPtr->count = found;

      // ++uRPtr;
    }else{

      uRPtr->count  = 0;  //needed if to few elements for acc >= ideal
      hisT   acc    = 0;
      hisT * hisIt  = histo.lbegin();

      // actual calculation of distribution
      for( T i = 0; i < histo.lsize()-1 ; ++i ) {
        acc += *hisIt++;

        if( acc >= ideal ){
          uRPtr->end       = i    ;
          uRPtr->count     = acc  ;
          acc              = 0    ;
          (++uRPtr)->begin = i + 1;
        }
      }
      if( acc >= MIN_NUM_ELEM_PER_UNIT ){ uRPtr->count = acc; }
      else
      {
        if( uRPtr == distr ){ uRPtr->count += acc; }
        else{ (--uRPtr)->count += acc; }
      }

      uRPtr->end = MAX_KEY ;

    }

    // set the rest to zero
    // while( ++uRPtr < distr_end ) { uRPtr->begin = 0; uRPtr->end = 0; uRPtr->count = 0; }
    if( ++uRPtr < distr_end){
      std::memset( uRPtr, 0, (distr_end - uRPtr) * sizeof(unitRange) );
    }

  } // end of unit 0 only part
  
 /* Distribute calculated distribution data to every unit.
  * Distribution data consists of unitRange structs.
  * With this every unit knows which unit is in charge for
  * sorting a specific value range.
  * Every unit is then creating buckets to send the values to
  * the responsible unit.
  */
  dart_ret_t ret = dart_bcast(
                      static_cast<void*>( distr )        ,  // buf
                      nUnits * sizeof(unitRange)         ,  // nelem
                      DART_TYPE_BYTE                     ,  // dtype
                      team_unit_t(0)                     ,  // root/source
                      team.dart_id( )                    ); // team

  if( DART_OK != ret ) cout << "An Error while BCAST has occured!" << endl;


 /* create dash::Array which will hold the data to be sorted.
  * therefore the CSRPattern of DASH will be used.
  * this allows different local sizes!
  */
  vector<extent_t> local_sizes;
  for( unitRange * uRPtr = distr; uRPtr < distr_end; ++uRPtr )
  {
      local_sizes.push_back( uRPtr->count );
  }


 /* Each unit creates buckets that are send to the responsible unit later.
  * Note: every unit will get data from other units (e.g. if distribution specifies range 0-0)
  * This is a array of vectors like "vector<Point> buckets[ involvedUnits ];"
  * But on the Heap!
  */
  const uint involvedUnits = local_sizes.size( );

  vector<Point> ** buckets = static_cast<vector<Point> **>(  std::malloc( involvedUnits * sizeof(vector<Point>*) )  );
  vector<Point> ** const buckets_end = buckets + involvedUnits;

  // create a vector for each bucket pointer
  for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++bucket ){ *bucket = new vector<Point>; }


  // iterate over pointsLocal and fill the buckets
  for( Point const * lclPt = pointsLocal.data(); lclPt < pointsLocal.data() + pointsLocal.size(); ++lclPt )
  {
    unitRange * uRPtr = distr;
    for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++uRPtr, ++bucket )
    {
      if( lclPt->value <= uRPtr->end )  // if value is in Range for this unit
      {
        (*bucket)->push_back(*lclPt);
        break;
      }
    }
  }

  // pointsLocal no longer needed
  pointsLocal.clear();
  pointsLocal.shrink_to_fit();

  pattern_t pattern( local_sizes );
  Array<Point, size_t, pattern_t> toSort( pattern );


  // not necessary anymore
  // std::memset(toSort.lbegin(), 0, toSort.lsize() * sizeof(Point));


 /* In comTable will be the information which unit has how many elements of type Point
  * for other units and itself.
  * The table can be thought of like "comTable[nUnits][local_sizes.size()]"
  * First every unit saves the information in row[0]
  */
  // these uint datatypes should be changed to long, but memory consumption with hundreds of untis
  // rises a lot and therefore this datatype has to suffice. otherwise memory consumption would be
  // that great that not enough elemts could be hold to need a long here
  uint * comTable = static_cast<uint*>(  std::malloc( sizeof(uint) * nUnits * involvedUnits )  );
  uint * const comTable_end = comTable + ( nUnits * involvedUnits );

  uint * thisUnitsRow_begin     = new uint[involvedUnits];
  uint * const thisUnitsRow_end = thisUnitsRow_begin + involvedUnits;

  uint * poiCount = thisUnitsRow_begin;
  for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++bucket, ++poiCount )
  {
    *poiCount = (*bucket)->size();
  }

  ret = dart_allgather(
           thisUnitsRow_begin,  // sendbuf
           comTable          ,  // recvbuf
           involvedUnits     ,  // nelem
           DART_TYPE_UINT    ,  // dtype
           team.dart_id( )   ); // team


  if( DART_OK != ret ) cout << "An Error while allgather has occured!" << endl;


 /* only involvedUnits can have data for memcpy!
  * and the units must have data for themselves as well
  * logical lookup -> comTable[myid][myid] > 0
  */
  uint forMySelf = comTable[myid * involvedUnits + myid];

  if( myid < involvedUnits && forMySelf > 0)
  {
    uint * forThisUnit = comTable + myid;
    hisT   lclOffset   = 0;


    for( uint ID = 0; ID < myid; ++ID)
    {
      // lclOffset += comTable[ID * involvedUnits + myid];  //replaced by pointer logic
      lclOffset   += *(forThisUnit);
      forThisUnit += involvedUnits ;
    }

    Point * lclDest = toSort.lbegin() + lclOffset;

    std::memcpy( lclDest, buckets[myid]->data(), sizeof(Point) * forMySelf );
  }


 /* copy data for other units in a circle.
  * every units sends its data to it's "right" neighbour
  */
  uint myRowOffset = myid * involvedUnits;

  auto baseIterator = toSort.begin( );
  // auto globDest = baseIterator +10;
  uint nextID = (myid + 1) % involvedUnits;


  for( uint counter = 0; counter < involvedUnits; ++counter, nextID = (nextID + 1) % involvedUnits )
  {
    // if this unit has data for the next unit and is a remote unit (local unit was handled before)
    if( comTable[ myRowOffset + nextID ] && nextID != myid  )
    {
          // uint   offsetToUnit  = 0;
          hisT   offsetOnUnit  = 0;
          uint * forNextUnit   = comTable + nextID;
      // extent_t * unitsLclSize  = local_sizes.data();

      // calculate offsetToUnit // not needed anymore since using .pattern.global(nextID, offsetOnUnit);
      // for( int ID = 0; ID < nextID; ++ID, ++unitsLclSize ){ offsetToUnit += *unitsLclSize; }
      // offsetOnUnit += comTable[ID * involvedUnits + nextID];  //replaced by pointer logic
      for( uint ID = 0; ID < myid; ++ID, forNextUnit += involvedUnits ){ offsetOnUnit += *(forNextUnit); }


      // global iterator/pointer to destination:
      // auto globDest = baseIterator + offsetToUnit + offsetOnUnit;
      auto globDest = baseIterator + toSort.pattern().global( team_unit_t(nextID), offsetOnUnit );

      //dash::copy / MPI_Put to target unit
      dash::copy( buckets[nextID]->data(), buckets[nextID]->data() + buckets[nextID]->size(), globDest );
    }
  }

  for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++bucket ){  delete *bucket; /* = new vector<Point>; */}

  // wait before sorting to finish all puts from dash::copy
  dash::barrier();


  std::stable_sort( toSort.lbegin( ), toSort.lend( ) );

  // }
 /* calculate local sizes for result array
  * local_sizes will be recycled therefore
  */

  local_sizes.clear();
  size_t chunk = toSort.size( ) / nelts;

  if( 0 == chunk ){

    local_sizes.push_back( nelts );
    for( uint id = 1; id < nUnits; ++id ) { local_sizes.push_back( 0 ); }

  }else
  {
    uint div        = 0     ;
    uint countSoFar = 0     ;
    uint lastIX     = 0     ;
    uint resulCount = 0     ;
    uint elemOnUnit = 0     ;
    uint remain     = nelts ;


    for( unitRange * uRPtr = distr; uRPtr < distr_end; ++uRPtr )
    {
        countSoFar += uRPtr->count;
        lastIX      = countSoFar - 1 ;
        div         = lastIX / chunk ;

        elemOnUnit = div + 1 - resulCount ;
        resulCount = div + 1              ;

        if( remain < elemOnUnit )
        {
          elemOnUnit = remain;
              remain = 0     ;
        }else
        {
          remain -= elemOnUnit;
        }

      local_sizes.push_back( elemOnUnit );

    }
  }

  // allocate result dash::Array with CSR Pattern
  pattern_t pattern_result  ( local_sizes    );
            result.allocate ( pattern_result );
  hisT pos = 0;

 /* generate local result.
  * start by calculating local start index.
  * using global indices therefore.
  */
  if( 0 < result.lsize( ) )
  {
   /* result.pattern().global(0) -> correspond to: "how much result elements have the units before this unit"
    * toSort.pattern().global(0) -> correspond to: "how much sorted elements have the units before this unit"
    */
         pos    = result.pattern().global(0) * chunk - toSort.pattern().global(0);
    Point * src = toSort.lbegin() + pos;

    for( value * dst = result.lbegin(); dst < result.lend(); ++dst, src += chunk )
    {
      dst->col = src->col;
      dst->row = src->row;
    }
  }

  dash::barrier();

  //no frees and deletes for potential speed up through cumulative memory freeing
}