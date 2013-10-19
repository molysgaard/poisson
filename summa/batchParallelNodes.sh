#!/bin/bash

# runs a batch of the paralell implementation for different node
# sets. Starts with 1 node, 1 mpi process, 2 cores and finishes at 6 nodes
# 6 mpi processes, 12 cores. 10 runs of each instance.
n_first=1
n_final=6
n_step=1

# run k repetitions of each job

k_first=1
k_final=10
k_step=1
 
for ((n=$n_first; n<=$n_final; n+=$n_step));
do
  for ((k=$k_first; k<=$k_final; k+=$k_step));
  do
    cat << EOF | qsub
#!/bin/bash -l
#PBS -N nodes
#PBS -l select=$n:ncpus=$((2*n)):mpiprocs=$n:mem=$((20/$n))G:cputype=E5-2670
#PBS -j oe
#PBS -l walltime=02:00:00

module load intel
module load mpt

cd $PBS_O_WORKDIR
cd summa

mpiexec_mpt ./summa_icc $n $((10000/$n)) 2
 
EOF

  done
done
