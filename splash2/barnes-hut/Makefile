TARGET=BARNES-dash.bin

DASHROOT?=$(HOME)/opt/dash-0.3.0
DASHCXX=$(DASHROOT)/bin/dash-mpic++ -dash:verbose

#CXXFLAGS= -ggdb3 -O0 -Wall -DENABLE_LOGGING -DENABLE_ASSERTIONS
CXXFLAGS= -O3 -Wall

.phony: all
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
