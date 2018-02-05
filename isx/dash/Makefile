
DASHROOT=$(HOME)/opt/dash-0.3.0/
DASHROOT=$(HOME)/opt-cray/dash-cce-debug/

#CC = mpicc
CC=cc
#CXX= $(DASHROOT)/bin/dash-mpic++ -dash:verbose
#CXX=mpicxx
CXX=CC
LD = $(CXX)
DEBUGFLAGS = -g -p -O0 -DDEBUG
OPTFLAGS = -O3 -DNDEBUG -mavx
CFLAGS += -Wall -std=c99 #$(OPTFLAGS)
CXXFLAGS += -Wall -std=c++11 -O3 -I$(DASHROOT)/include #-DUSE_ATOMICS
LDLIBS += -lrt -lm -L$(DASHROOT)/lib -ldash-mpi -ldart-mpi -ldart-base #-lpapi -lhwloc -lnuma
LDFLAGS += 
### 
# If you change these values, you must also change them in params.h
STRONG = 1
WEAK = 2
WEAKISO = 3
###

.PHONY: all clean debug

SRCS := $(wildcard *.cc)
OBJDIR := obj
BINDIR := bin
STRONG_EXE = isx.strong
WEAK_EXE = isx.weak 
WEAKISO_EXE = isx.weak_iso
STRONG_OBJS := $(SRCS:.cc=.o_s)
STRONG_OBJS := $(addprefix $(OBJDIR)/,$(STRONG_OBJS))
WEAK_OBJS := $(SRCS:.cc=.o_w)
WEAK_OBJS := $(addprefix $(OBJDIR)/,$(WEAK_OBJS))
WEAKISO_OBJS := $(SRCS:.cc=.o_wi)
WEAKISO_OBJS := $(addprefix $(OBJDIR)/,$(WEAKISO_OBJS))

all: $(STRONG_EXE) $(WEAK_EXE) $(WEAKISO_EXE)

optimized: CFLAGS += $(OPTFLAGS)
optimized: all

debug: CFLAGS += $(DEBUGFLAGS)
debug: all

$(STRONG_EXE):$(STRONG_OBJS)
	@mkdir -p $(BINDIR)
	$(LD) $(STRONG_OBJS) -o $(BINDIR)/$(STRONG_EXE) $(LDLIBS) $(LDFLAGS)

$(WEAK_EXE):$(WEAK_OBJS)
	@mkdir -p $(BINDIR)
	$(LD) $(WEAK_OBJS) -o $(BINDIR)/$(WEAK_EXE) $(LDLIBS) $(LDFLAGS)

$(WEAKISO_EXE):$(WEAKISO_OBJS)
	@mkdir -p $(BINDIR)
	$(LD) $(WEAKISO_OBJS) -o $(BINDIR)/$(WEAKISO_EXE) $(LDLIBS) $(LDFLAGS)

$(OBJDIR)/%.o_s : %.cc
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -D SCALING_OPTION=$(STRONG) -c $< -o $@

$(OBJDIR)/%.o_w: %.cc
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -D SCALING_OPTION=$(WEAK) -c $< -o $@

$(OBJDIR)/%.o_wi: %.cc
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -D SCALING_OPTION=$(WEAKISO) -c $< -o $@

clean:
	@- $(RM) $(BINDIR)/$(STRONG_EXE)
	@- $(RM) $(STRONG_OBJS)
	@- $(RM) $(BINDIR)/$(WEAK_EXE)
	@- $(RM) $(WEAK_OBJS)
	@- $(RM) $(BINDIR)/$(WEAKISO_EXE)
	@- $(RM) $(WEAKISO_OBJS)

depend:
	    makedepend -- $(CFLAGS) -- $(SRCS)

# Don't place anything below this line, since
# # the make depend program will overwrite it
# # DO NOT DELETE THIS LINE -- make depend depends on it.
