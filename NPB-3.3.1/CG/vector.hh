#ifndef VECTOR_HH_
#define VECTOR_HH_

#include <libdash.h>

template<class TValue>
class Vector {

  dash::Array<TValue> _internal;

public:
  Vector(int size) : _internal(size) {
  }
};

#endif // VECTOR_HH_
