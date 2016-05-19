#ifndef LULESH_UTIL_H_INCLUDED
#define LULESH_UTIL_H_INCLUDED

#include <iostream>
#include "lulesh.h"
#include "lulesh-dash.h"

template<typename T>
void peek(T *ptr, size_t nval, std::ostream& os=std::cout)
{
  size_t i;
  for( i=0; i<nval; ++i ) {
    os << (*ptr++) << " ";
  }
  os << std::endl;
}

template<typename T>
T chksum(T *ptr, size_t nval)
{
  T res = 0.0;
  size_t i;
  for( i=0; i<nval; ++i ) {
    res += (i+1) * (*ptr++);
  }
  return res;
}

void print_config(Domain& dom,
		  std::ostream& os=std::cout);

#endif // LULESH_UTIL_H_INCLUDED
