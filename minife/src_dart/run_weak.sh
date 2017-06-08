#!/bin/bash

#PBS -lnodes=344:ppn=24
#PBS -lwalltime=2:00:00
#PBS -n 
#PBS -N minife_rma
#PBS -m ae
#PBS -j oe

cd  ~/src/dash/dash-apps/minife/src_dart

module swap PrgEnv-cray PrgEnv-gnu

num_nodes=${PBS_NUM_NODES}

# run on 1-8k processor (full nodes)
for f in 1 2 4 8 ; do
	echo "[$(date +%H:%M:%S)] Preparing for f=$f"
	# 1k ~ 1032 procs ~ 43*24
	num_nodes=$((43 * f))
	num_procs=$((num_nodes*24))
	x=1000
	y=$x
	z=$((100*f))
	echo "[$(date +%H:%M:%S)] Running with $num_procs processors ($x $y $z)"

	# run without DMAPP first
	unset MPICH_RMA_OVER_DMAPP

	# DART
	echo "[$(date +%H:%M:%S)] Running DART dynamic (nodmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.dynamic.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1  &> miniFE.weak.dynamic.nodmapp.$num_procs.log
	echo "[$(date +%H:%M:%S)] Running DART allocate (nodmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.allocate.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1 &> miniFE.weak.allocate.nodmapp.$num_procs.log

	# DASH
	echo "[$(date +%H:%M:%S)] Running DASH (nodmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.dash.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1  &> miniFE.weak.dash.nodmapp.$num_procs.log
	
	# MPI
	echo "[$(date +%H:%M:%S)] Running MPI-RMA (nodmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.mpirma.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1  &> miniFE.weak.mpirma.nodmapp.$num_procs.log
	echo "[$(date +%H:%M:%S)] Running MPI (nodmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.ref.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1 &> miniFE.weak.ref.nodmapp.$num_procs.log

	# enable DMAPP
	export MPICH_RMA_OVER_DMAPP=1

	# DART 
	echo "[$(date +%H:%M:%S)] Running DART allocate (dmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.dynamic.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1 &> miniFE.weak.dynamic.dmapp.$num_procs.log
	echo "[$(date +%H:%M:%S)] Running DART dynamic (dmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.allocate.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1 &> miniFE.weak.allocate.dmapp.$num_procs.log

	# DASH
	echo "[$(date +%H:%M:%S)] Running DASH (dmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.dash.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1 &> miniFE.weak.dash.dmapp.$num_procs.log

	# MPI
	echo "[$(date +%H:%M:%S)] Running MPI-RMA (dmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.mpirma.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1  &> miniFE.weak.mpirma.dmapp.$num_procs.log
	echo "[$(date +%H:%M:%S)] Running MPI (dmapp)"
	time aprun -n $num_procs -N 24 ./miniFE.ref.x -nx $x -ny $y -nz $z -mv_overlap_comm_comp 1 -stdout 1 &> miniFE.weak.ref.dmapp.$num_procs.log
	echo "[$(date +%H:%M:%S)] Done."
done
