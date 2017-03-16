/* #ifndef DASH_ENABLE_LOGGING
#define DASH_ENABLE_LOGGING
#endif */

#include <libdash.h>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

using uint  = unsigned int;

inline void product(         
                     const uint   nelts ,
           vector<double> const & vec,
  const dash::NArray<double, 2> & matIn,
            dash::Array<double> & result,
                      const int   myid
    );

template< typename T >
void print2d( T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>( mat(i,j) )<< " ";
    }
    cout << endl;
  }
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

int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  int myid = static_cast<int>( dash::Team::GlobalUnitID( ).id );
  
  
  if( argc != 2 ){
    if( 0 == myid ){ cout << "1 Parameter expected!\n"
                          << "Usage: cowichan_product nElements\n"
                          << "Then enter the matrix and the vector array." << endl;
    }
    dash::finalize( );
    return 0;
  }
  
  uint nelts = static_cast<uint>( atoi( argv[1] ) );
  
  dash::NArray<double, 2> matIn(nelts, nelts);
  dash::Array <double>    result(nelts);
  vector <double>         vec(nelts);
  
  //read input on unit 0
  if( 0 == myid ) {
    readMatrix(matIn);
    readVec(vec);
    cout << endl;
  }
  
  //broadcast input from unit 0
  broadcastInputToUnits(vec);

  product(nelts, vec, matIn, result, myid);
  
  dash::barrier();
  
  if(0 == myid){
    cout << nelts << "\n";
    for(auto i : result) {
      cout << static_cast<double>(i) << " ";
    } cout << endl;
  }
  
  dash::finalize( );
}


inline void product(         
                     const uint   nelts ,
           vector<double> const & vec,
  const dash::NArray<double, 2> & matIn,
            dash::Array<double> & result,
                      const int   myid
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
