include make.defs

.phony: all
all: $(TARGET)

OBJS = code.o getparam.o load.o grav.o util.o

.cc.o:
	$(CXX) $(INC) $(CXXFLAGS) -c $*.cc

.c.o:
	$(CC) $(INC) -c $(CFLAGS) $*.c

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(CXXFLAGS) $(LDFLAGS) $(LIB) -o $(TARGET)

.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean

clean:
	rm -f $(TARGET) *.o
