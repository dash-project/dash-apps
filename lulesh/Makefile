
include make.defs

#CFLAGS += -DUSE_MPI
CFLAGS += -DUSE_DASH


lulesh: lulesh.o lulesh-opts.o lulesh-dash.o lulesh-util.o	\
	lulesh-calc.o lulesh-dash-regions.o lulesh-comm-mpi.o	\
	lulesh-comm-mpi-sendrecv.o lulesh-comm-dash.o		\
	lulesh-comm-dash-onesided.o
	$(CXX) $(LFLAGS) -o $@ $+ $(LIBDASH) $(LIBDART) $(LIBDART_EXTRA)

lulesh-dash : lulesh.o lulesh-opts.o lulesh-dash.o lulesh-util.o	\
	lulesh-calc.o lulesh-comm-mpi.o lulesh-comm-mpi-sendrecv.o	\
	lulesh-dash-regions.o
	$(CXX) $(LFLAGS) -o $@ $+ $(LIBDASH) $(LIBDART) $(LIBDART_EXTRA)

%.o 	: %.cc %.h
	$(CXX) -c $(CFLAGS) -I$(DARTIF_INC) -I$(DART_INC) -I$(DASH_INC) $<


printenv :
	@echo "CXX           = $(CXX)"
	@echo "CFLAGS        = $(CFLAGS)"
	@echo "LFLAGS        = $(LFLAGS)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "LIBDASH       = $(LIBDASH)"
	@echo "LIBDART       = $(LIBDART)"
	@echo "LIBDART_EXTRA = $(LIBDART_EXTRA)"
	@echo "DARTIF_INC    = $(DARTIF_INC)"
	@echo "DART_INC      = $(DART_INC)"
	@echo "DASH_INC      = $(DASH_INC)"


clean :
	rm -f *~
	rm -f *.o
	rm -f *.gch
	rm -f ./lulesh

