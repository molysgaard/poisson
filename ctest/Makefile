ctest : ctest.c
	mpicc ctest.c -o ctest

run : ctest
	mpirun -n 2 ctest
