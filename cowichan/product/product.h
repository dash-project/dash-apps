using std::cout;
using std::endl;
using std::vector;

using dash::NArray;
using dash::Array;
using dash::team_unit_t;
using dash::Team;
using dash::barrier;

inline void PrintOutput( Array <double> const & result )
{
  if(0 == myid){
    cout << nelts << "\n";
    cout << std::showpoint << std::fixed << std::setprecision(4);
    
    for(auto i : result) {
      cout << static_cast<double>(i) << " ";
    } cout << endl;
  }
}


inline void Product(
  vector < double    > const & vec    ,
  NArray < double, 2 > const & matIn  ,
  Array  < double    >       & result )
{
   uint lclRows = matIn.pattern().local_extents()[0];
   double         sum  ;
   double const * mPtr ;
   double const * vPtr ;
   double       * res = result.lbegin();

   for(uint i = 0; i < lclRows; ++i){
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
