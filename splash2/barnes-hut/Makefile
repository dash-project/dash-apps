TARGET=BARNES-dash.bin

DASHROOT?=$(HOME)/opt/dash-0.3.0
DASHCXX=$(DASHROOT)/bin/dash-mpic++ -dash:verbose

#CXXFLAGS= -ggdb3 -O0 -Wall -DENABLE_LOGGING -DENABLE_ASSERTIONS
CXXFLAGS= -O3 -Wall

DEBUG_RUN=1

.phony: all run
all: $(TARGET)

OBJS = code.o getparam.o load.o grav.o util.o

.cc.o:
	$(DASHCXX) $(CXXFLAGS) -c $*.cc

.c.o:
	$(CC) -c $(CFLAGS) $*.c

$(TARGET): $(OBJS)
	$(DASHCXX) $(OBJS) $(CXXFLAGS) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o

run:
	$(eval $@_TMP := $(shell mktemp /tmp/nbody.XXXXXXXX))
	mpirun -n 4 ./$(TARGET) 64 $(DEBUG_RUN) > $($@_TMP)
	./verify/verify.sh $($@_TMP) verify/barnes.64.gold

