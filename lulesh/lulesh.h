#ifndef LULESH_H_INCLUDED
#define LULESH_H_INCLUDED

#include <iostream>

// Precision specification
typedef float        real4;
typedef double       real8;
typedef long double  real10;  // 10 bytes on x86

typedef int    Index_t;  // array subscript and loop index
typedef real8  Real_t;   // floating point representation
typedef int    Int_t;    // integer representation

template<typename T>
void peek(T *ptr, int len, std::ostream& os=std::cout);

#endif // LULESH_H_INCLUDED

