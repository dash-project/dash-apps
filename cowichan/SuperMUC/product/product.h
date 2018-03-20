using std::cout;
using std::endl;
using std::vector;

using dash::NArray;
using dash::Array;
using dash::team_unit_t;
using dash::Team;
using dash::barrier;

inline void PrintOutput( Array <double> const & result, uint const nelts)
{
  if(0 == myid){
    cout << nelts << "\n";
    cout << std::showpoint << std::fixed << std::setprecision(4);

    for(auto i : result) {
      cout << static_cast<double>(i) << " ";
    } cout << endl;
  }
}


inline void BroadcastOuterVecToUnits( vector <double> & vec )
{
  team_unit_t TeamUnit0ID = Team::All().myid( );
  TeamUnit0ID.id = 0;
  dart_ret_t ret = dart_bcast(
                    static_cast<void*>(vec.data( )),  // buf
                    vec.size( )                    ,  // nelts
                    DART_TYPE_DOUBLE               ,  // dtype
                    TeamUnit0ID                    ,  // root
                    Team::All().dart_id( )         ); // team

  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl;
}


inline void Product(
  vector < double    > const & vec    ,
  NArray < double, 2 > const & matIn  ,
  Array  < double    >       & result ,
                 uint  const   nelts  )
{
   uint lclRows = matIn.pattern().local_extents()[0];
   double         sum  ;
   double const * mPtr ;
   double const * vPtr ;
   double       * res = result.lbegin();

   for(uint i = 0; i < lclRows; ++i)
   {
     sum  = 0;
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

   barrier();
}
