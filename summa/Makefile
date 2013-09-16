summa: summa.c
	mpicc -O3 summa.c -lblas -llapack -o summa

run: summa
	mpirun -n 4 summa

clean:
	rm summa
