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
