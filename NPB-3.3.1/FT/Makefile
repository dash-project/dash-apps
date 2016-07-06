SHELL=/bin/sh
BENCHMARK=ft
BENCHMARKU=FT

include ../config/make.def

include ../sys/make.common
#Override PROGRAM
FTPROGRAM  = $(BINDIR)/$(BENCHMARK).$(CLASS).x

OBJS = ft.o \
	${COMMON}/c_print_results.o ${COMMON}/c_timers.o ${COMMON}/c_randdp.o

${PROGRAM}: config bindir ${OBJS}
	${CXXLINK}  -o ${FTPROGRAM} ${OBJS} ${CMPI_LIB} ${CLINKFLAGS}

bindir:
	test -d $(BINDIR) || mkdir $(BINDIR)

%.o : %.cpp
	${CXXCOMPILE} $<

.c.o:
	${CCOMPILE} $<

ft.o:             ft.cpp  npbparams.h

clean:
	- rm -f *.o *~ mputil*
	- rm -f ft npbparams.h core
