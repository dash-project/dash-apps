/* #ifndef DASH_ENABLE_LOGGING
#define DASH_ENABLE_LOGGING
#endif */

#include <libdash.h>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

using uint  = unsigned int;

uint nelts;
static int myid;


template< typename T >
inline void PrintOutput( T const & result ) {
  if(0 == myid){
    cout << nelts << "\n";
    cout << std::showpoint << std::fixed << std::setprecision(4);
    for(auto i : result) {
      cout << static_cast<double>(i) << " ";
    } cout << endl;
  }
}


template<typename T, typename X>
inline void ReadMatrixAndVector(T& matIn, X& vec){
  if( 0 == myid ) {
    //Read Matrix
    double tmp;
    for ( auto i : matIn ){
      cin >> tmp;
      i = tmp;
    }
    
    //Read Vector
    for (int i = 0; i < vec.size(); i++) {
      cin >> vec[i];
    }
  }
}


inline void ReadNelts( ){
  
  dash::Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


inline void BroadcastInputToUnits(vector<double> & vec) {
  dash::team_unit_t TeamUnit0ID = dash::Team::All().myid( );
  TeamUnit0ID.id = 0;
  dart_ret_t ret = dart_bcast(
                      static_cast<void*>(vec.data( )),  // buf 
                      vec.size( )                    ,  // nelts
                      DART_TYPE_DOUBLE               ,  // dtype
                      TeamUnit0ID                    ,  // root
                      dash::Team::All().dart_id( )      // team
                   );
  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl; 
}


template<typename T, typename X, typename Y>
inline void Product(
     X const & vec,
     T const & matIn,
           Y & result
    ){
   uint lclRows = matIn.pattern().local_extents()[0];
   double * res = result.lbegin();
   double   sum;
   double const * mPtr;
   double const * vPtr;

   for(uint i = 0; i < lclRows; ++i){
     sum = 0;
     mPtr = matIn.local.row(i).lbegin();
     vPtr = vec.data();
     
     /* first loop iteration is done here.
      * so in loop can the prefix operator be used
      */
      sum += *mPtr * *vPtr;
     
     /* one less loop iteration because of the line before
      * -> j != 0 but j = 1
      */
     for(uint j = 1; j < nelts ; ++j){
       sum += *++mPtr * *++vPtr;
     }
     *(res++) = sum;
   }
   
   dash::barrier();
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  myid = dash::myid( );
  
  ReadNelts( );
  
  dash::NArray<double, 2> matIn(nelts, nelts);
  dash::Array <double>    result(nelts);
  vector <double>         vec(nelts);
 
  //read input on unit 0 and broadcast it
  ReadMatrixAndVector(matIn, vec);
  BroadcastInputToUnits(vec);

  Product(vec, matIn, result);
  PrintOutput( result );
  
  dash::finalize( );
}