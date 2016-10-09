#ifndef VECTOR_HH_
#define VECTOR_HH_

#include <libdash.h>

template<class TValue>
class Vector {

  dash::Array<TValue> & _internal;
  
public:
  Vector(dash::Array<TValue> & internal) : _internal(internal) { }

  TValue Norm() {

    dash::Array<TValue> l_results(dash::size());
    l_results.local[0] =
      std::accumulate(
		       _internal.lbegin(),
		       _internal.lend(),
		       0,
		       [](TValue x1, TValue x2){
			 return x1 + x2*x2;
		       });


    return  dash::accumulate(l_results.begin(),
			     l_results.end(),
			     0, dash::plus<TValue>());
  }

};

#endif // VECTOR_HH_
