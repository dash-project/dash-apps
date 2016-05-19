#ifndef LULESH_UTIL_H_INCLUDED
#define LULESH_UTIL_H_INCLUDED

#include <iostream>
#include "lulesh.h"

template<typename T>
void peek(T *ptr, size_t nval, std::ostream& os=std::cout);

template<typename T>
T chksum(T *ptr, size_t nval);

void print_config(const Domain& dom);

#endif // LULESH_UTIL_H_INCLUDED
