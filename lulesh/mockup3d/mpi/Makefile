
CXX = mpicxx -Wall

all : lulesh-mock.o domain.o comm.o
	$(CXX) -o lulesh-mock $^

%.o : %.cc
	$(CXX) -c -o $@ $+

clean :
	rm -f ./*.o
	rm -f *~
	rm -f ./lulesh-mock
