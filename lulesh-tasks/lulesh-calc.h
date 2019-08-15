#ifndef LULESH_CALC_H_INCLUDED
#define LULESH_CALC_H_INCLUDED

#include "lulesh.h"
#include "lulesh-dash.h"

//
// LULESH "calculation" routines, largely unchanged from their
// original version in the MPI version of LULESH
//
void CalcTimeConstraintsForElems(Domain& domain);


void CalcKinematicsForElems(Domain &domain,
			    Real_t deltaTime, Index_t numElem);

Real_t CalcElemVolume(const Real_t x[8],
		      const Real_t y[8],
		      const Real_t z[8]);

void CalcVolumeForceForElems(Domain& domain);

void CalcMonotonicQForElems(Domain& domain);

void CalcMonotonicQGradientsForElems(Domain& domain);

void EvalEOSForElems(Domain& domain, Real_t *vnewc,
                     Int_t numElemReg,
		     Index_t *regElemList, Int_t rep);

#endif /* LULESH_CALC_H_INCLUDED */
