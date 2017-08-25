all: main
	./$< < $<.in > main.out
	diff main.out main.gold
	./$< < $<2.in > main2.out
	diff main2.out main2.gold
	./$< < $<3.in > main3.out
	diff main3.out main3.gold
	./$< < $<4.in > main4.out
	diff main4.out main4.gold
	./$< < $<5.in > main5.out
	diff main5.out main5.gold
	./$< < $<6.in > main6.out
	diff main6.out main6.gold
	./main --is_bench < main.in 
main: main.cpp
	icpc ${CFLAGS} -O2 -o main main.cpp
