SHELL=/bin/sh
BENCHMARK=dt
BENCHMARKU=DT

include ../config/make.def

include ../sys/make.common
#Override PROGRAM
DTPROGRAM  = $(BINDIR)/$(BENCHMARK).$(CLASS).x

OBJS = dt.o DGraph.o \
	${COMMON}/c_print_results.o ${COMMON}/c_timers.o ${COMMON}/c_randdp.o


${PROGRAM}: config bindir ${OBJS}
	${CXXLINK}  -o ${DTPROGRAM} ${OBJS} ${CMPI_LIB} ${CLINKFLAGS}

bindir:
	test -d $(BINDIR) || mkdir $(BINDIR)

%.o : %.cpp
	${CXXCOMPILE} $<

.c.o:
	${CCOMPILE} $<

dt.o:             dt.cpp  npbparams.h
DGraph.o:	DGraph.cpp DGraph.h

clean:
	- rm -f *.o *~ mputil*
	- rm -f dt npbparams.h core
