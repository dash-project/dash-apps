#ifndef LULESH_H_INCLUDED
#define LULESH_H_INCLUDED

#include <cmath>
#include <iostream>

// Precision specification
typedef float        real4;
typedef double       real8;
typedef long double  real10;  // 10 bytes on x86

typedef int    Index_t;  // array subscript and loop index
typedef real8  Real_t;   // floating point representation
typedef int    Int_t;    // integer representation

enum { VolumeError = -1, QStopError = -2 } ;

inline real4  SQRT(real4  arg) { return sqrtf(arg) ; }
inline real8  SQRT(real8  arg) { return sqrt(arg) ; }
inline real10 SQRT(real10 arg) { return sqrtl(arg) ; }

inline real4  CBRT(real4  arg) { return cbrtf(arg) ; }
inline real8  CBRT(real8  arg) { return cbrt(arg) ; }
inline real10 CBRT(real10 arg) { return cbrtl(arg) ; }

inline real4  FABS(real4  arg) { return fabsf(arg) ; }
inline real8  FABS(real8  arg) { return fabs(arg) ; }
inline real10 FABS(real10 arg) { return fabsl(arg) ; }

template<typename T>
void peek(T *ptr, int len, std::ostream& os=std::cout);

#endif // LULESH_H_INCLUDED

