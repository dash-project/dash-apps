#ifndef LULESH_CALC_H_INCLUDED
#define LULESH_CALC_H_INCLUDED

#include "lulesh-dash.h"

//
// LULESH "calculation" routines, largely unchanged from their
// original version in the MPI version of LULESH
//
void CalcTimeConstraintsForElems(Domain& domain);


void CalcKinematicsForElems(Domain &domain, Real_t *vnew,
			    Real_t deltaTime, Index_t numElem);

Real_t CalcElemVolume(const Real_t x[8],
		      const Real_t y[8],
		      const Real_t z[8]);

void CalcVolumeForceForElems(Domain& domain);



#endif /* LULESH_CALC_H_INCLUDED */
