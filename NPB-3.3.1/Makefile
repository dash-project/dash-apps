SHELL=/bin/sh
CLASS=U
NPROCS=1
SUBTYPE=
VERSION=
SFILE=config/suite.def

default: header
	@ sys/print_instructions

DT: dt
dt: header
	cd DT; $(MAKE) CLASS=$(CLASS)

dtS: clean
	cd DT; $(MAKE) NPROCS=10 CLASS=S

dtW: clean
	cd DT; $(MAKE) NPROCS=10 CLASS=W

dtA: clean
	cd DT; $(MAKE) NPROCS=10 CLASS=A
# Awk script courtesy cmg@cray.com, modified by Haoqiang Jin
suite:
	@ awk -f sys/suite.awk SMAKE=$(MAKE) $(SFILE) | $(SHELL)


# It would be nice to make clean in each subdirectory (the targets
# are defined) but on a really clean system this will won't work
# because those makefiles need config/make.def
clean:
	- rm -f core 
	- rm -f *~ */core */*~ */*.o */npbparams.h */*.obj */*.exe
	- rm -f MPI_dummy/test MPI_dummy/libmpi.a
	- rm -f sys/setparams sys/makesuite sys/setparams.h
	- rm -f btio.*.out*

veryclean: clean
	#- rm -f config/make.def config/suite.def 
	- rm -f bin/sp.* bin/lu.* bin/mg.* bin/ft.* bin/bt.* bin/is.* 
	- rm -f bin/ep.* bin/cg.* bin/dt.*

header:
	@ sys/print_header

