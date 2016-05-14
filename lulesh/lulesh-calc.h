#ifndef LULESH_CALC_H_INCLUDED
#define LULESH_CALC_H_INCLUDED

static inline
Real_t CalcElemVolume( const Real_t x0, const Real_t x1,
		       const Real_t x2, const Real_t x3,
		       const Real_t x4, const Real_t x5,
		       const Real_t x6, const Real_t x7,
		       const Real_t y0, const Real_t y1,
		       const Real_t y2, const Real_t y3,
		       const Real_t y4, const Real_t y5,
		       const Real_t y6, const Real_t y7,
		       const Real_t z0, const Real_t z1,
		       const Real_t z2, const Real_t z3,
		       const Real_t z4, const Real_t z5,
		       const Real_t z6, const Real_t z7 );

Real_t CalcElemVolume( const Real_t x[8],
		       const Real_t y[8],
		       const Real_t z[8] );

#endif /* LULESH_CALC_H_INCLUDED */
