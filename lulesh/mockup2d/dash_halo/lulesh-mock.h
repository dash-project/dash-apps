#ifndef LULESH_MOCK_H_INCLUDED
#define LULESH_MOCK_H_INCLUDED

using Real_t = double;

struct TripleT {
  Real_t x;
  Real_t y;
  Real_t z;

  TripleT& operator += ( const TripleT other) {
    x += other.x;
    y += other.y;
    z += other.z;
  }

  TripleT operator + ( const TripleT other) {
    TripleT new_triple(*this);
    new_triple += other;

    return new_triple;
  }
};

typedef int     Index_t;

#endif /* LULESH_MOCK_H_INCLUDED */
