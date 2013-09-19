compile: summa sequential

summa: summa.c
	mpicc -O3 summa.c -lblas -llapack -lm -o summa

sequential: sequential.c
	gcc -O3 sequential.c -llapack -lblas -lm -lrt -o sequential


run: summa
	mpirun -n 4 summa

clean:
	rm summa
