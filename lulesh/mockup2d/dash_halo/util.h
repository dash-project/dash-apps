#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <iostream>
#include <iomanip>

template<typename T>
void print2d(const T& mat, std::ostream& os)
{
  for( int j = 0; j < mat.extent(1); j++ ) {
    for( int i = 0; i < mat.extent(0); i++ ) {
      os << std::setw(3) <<
	static_cast<typename T::value_type>(mat(i,j)) << " ";
    }
    os << std::endl;
  }
}

template<typename T>
void print2dTriple(const T& mat, std::ostream& os)
{
  for( int j = 0; j < mat.extent(1); j++ ) {
    for( int i = 0; i < mat.extent(0); i++ ) {
      const auto& triple = static_cast<typename T::value_type>(mat(i,j));
      os << std::setw(5)
         << triple.x << "," << std::setw(2)
         << triple.y << "," << std::setw(2)
         << triple.z << " ";
    }
    os << std::endl;
  }
}
#endif /* UTIL_H_INCLUDED */

