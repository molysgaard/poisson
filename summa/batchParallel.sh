#!/bin/bash

# runs a batch of the paralell implementation for different problem sizes
# with a fixed 5 node, 5 mpiprocs, 2 openmp threads
 
n_first=1000
n_final=10000
n_step=1000

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
#PBS -N probsize
#PBS -l select=5:ncpus=10:mpiprocs=5:mem=3G:cputype=E5-2670
#PBS -j oe
#PBS -l walltime=00:10:00

module load intel
module load mpt

cd $PBS_O_WORKDIR
cd summa

mpiexec_mpt ./summa_icc 5 $(($n/5)) 2
 
EOF

  done
done
