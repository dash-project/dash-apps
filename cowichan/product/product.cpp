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
void print2d( T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>( mat(i,j) )<< " ";
    }
    cout << endl;
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


inline void readMatrix(dash::NArray<double, 2> & matIn) {
    double tmp;
    for ( auto i : matIn ){
      cin >> tmp;
      i = tmp;
    }
}


 inline void readVec( vector<double> & vec ) {
  for (int i = 0; i < vec.size(); i++) {
    cin >> vec[i];
  }
}


inline void broadcastInputToUnits(vector<double> & vec) {
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


inline void product(
           vector<double> const & vec,
  dash::NArray<double, 2> const & matIn,
            dash::Array<double> & result
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
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  myid = dash::myid( );
  
  ReadNelts( );
  
  dash::NArray<double, 2> matIn(nelts, nelts);
  dash::Array <double>    result(nelts);
  vector <double>         vec(nelts);
  
  //read input on unit 0
  if( 0 == myid ) {
    readMatrix(matIn);
    readVec(vec);
  }
  
  //broadcast vector input from unit 0
  broadcastInputToUnits(vec);

  product(vec, matIn, result);
  
  dash::barrier();
  
  if(0 == myid){
    cout << nelts << "\n";
    cout << std::showpoint << std::fixed << std::setprecision(4);
    for(auto i : result) {
      cout << static_cast<double>(i) << " ";
    } cout << endl;
  }
  
  dash::finalize( );
}