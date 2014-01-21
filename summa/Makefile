compile: summa sequential

summa_icc: summa.c
	#icc -O2 -openmp -mkl=parallel -lmpi summa.c -o summa_icc
	icc -O2 summa.c -o summa_icc  -DMKL_ILP64 -I$(MKLROOT)/include -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -lmpi

sequential_mkl: sequential_mkl.c
	#icc -O2 sequential_mkl.c -o sequential_mkl -DMKL_ILP64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -lrt
	icc -O2 sequential_mkl.c -o sequential_mkl -DMKL_ILP64 -I$(MKLROOT)/include -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -lrt

summa: summa.c
	mpicc -O3 summa.c -fopenmp -lblas -lm -lrt -o summa

sequential: sequential.c
	gcc -O2 sequential.c -lblas -lm -lrt -o sequential

sequential_icc: sequential.c
	icc -O2 sequential.c -lm -lrt -o sequential_icc

clean:
	rm summa summa_icc sequential sequential_mkl sequential_icc u.m u-seq.m poissonPar.eps
