#!/bin/bash

# runs a batch of the serial implementation for different problem sizes
 
n_first=1000
n_final=10000
n_step=1000

# number of repetitions of each benchmark
k_first=1
k_final=10
k_step=1
 
for ((n=$n_first; n<=$n_final; n+=$n_step));
do
  for ((k=$k_first; k<=$k_final; k+=$k_step));
  do
    cat << EOF | qsub
#!/bin/bash -l
#PBS -N ser$n-$k
#PBS -l select=1:ncpus=1:mem=3G:cputype=E5-2670
#PBS -j oe
#PBS -l walltime=03:00:00

cd $PBS_O_WORKDIR
cd summa

./sequential_mkl $n
 
EOF
 
  done
done
