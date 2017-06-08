#!/bin/bash

#PBS -lnodes=344:ppn=24
#PBS -lwalltime=2:00:00
#PBS -n 
#PBS -N hpccg_rma
#PBS -m ae
#PBS -j oe

cd ~/src/dash/dash-apps/hpccg/src_dart

module swap PrgEnv-cray PrgEnv-gnu

num_nodes=${PBS_NUM_NODES}

# run on 1-8k processor (full nodes)
for f in 1 2 4 8 ; do
	# 1k ~ 1032 procs == 43*24
	num_nodes=$((43 * f))
	num_procs=$((num_nodes*24))
	# HPCCG takes size arguments per-process, scales automatically
	x=100
	y=$x
	z=$x
	echo "Running with $num_procs processors ($x $y $z)"

	# run without DMAPP first
	unset MPICH_RMA_OVER_DMAPP

	# DART
	time aprun -n $num_procs -N 24 ./HPCCG.dynamic.x $x $y $z  &> HPCCG.weak.dynamic.nodmapp.$num_procs.log
	time aprun -n $num_procs -N 24 ./HPCCG.allocate.x $x $y $z &> HPCCG.weak.allocate.nodmapp.$num_procs.log
	
	# DASH
	time aprun -n $num_procs -N 24 ./HPCCG.dash.x $x $y $z &> HPCCG.weak.dash.nodmapp.$num_procs.log

	# MPI
	time aprun -n $num_procs -N 24 ./HPCCG.mpirma.x $x $y $z &> HPCCG.weak.mpirma.nodmapp.$num_procs.log
	time aprun -n $num_procs -N 24 ./HPCCG.ref.x $x $y $z    &> HPCCG.weak.ref.nodmapp.$num_procs.log

	# enable DMAPP
	export MPICH_RMA_OVER_DMAPP=1

	# DART
	time aprun -n $num_procs -N 24 ./HPCCG.dynamic.x $x $y $z  &> HPCCG.weak.dynamic.dmapp.$num_procs.log
	time aprun -n $num_procs -N 24 ./HPCCG.allocate.x $x $y $z &> HPCCG.weak.allocate.dmapp.$num_procs.log
	
	# DASH
	time aprun -n $num_procs -N 24 ./HPCCG.dash.x $x $y $z &> HPCCG.weak.dash.dmapp.$num_procs.log

	# MPI
	time aprun -n $num_procs -N 24 ./HPCCG.mpirma.x $x $y $z  &> HPCCG.weak.mpirma.dmapp.$num_procs.log
	time aprun -n $num_procs -N 24 ./HPCCG.ref.x $x $y $z     &> HPCCG.weak.ref.dmapp.$num_procs.log
done
