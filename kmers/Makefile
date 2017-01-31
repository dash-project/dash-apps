include make.defs

KMER_LENGTH 		= 19
KMER_PACKED_LENGTH 	= $(shell echo $$((($(KMER_LENGTH)+3)/4)))

CFLAGS = -Iinclude

.phony: all
all: serial dash

serial: CFLAGS += -O2 -std=c99 -g -D_POSIX_C_SOURCE=1
serial: serial.c
	$(CC) $(CFLAGS) -o $@ $< -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH) $(LIBS)

dash:  dash.cpp
	$(CXX) -c $(INC) $?
	$(CXX) -o $@ $@.o $(LIB)

.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f serial dash *.o
