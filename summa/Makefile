summa: summa.c
	mpicc -O2 summa.c -lblas -llapack -o summa

run: summa
	mpirun -n 4 summa

clean:
	rm summa
