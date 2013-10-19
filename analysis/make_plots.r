#!/usr/bin/Rscript

library(gdata)

#pdf('serial_walltime.pdf')
#serData = read.csv('ser-1k-10k-20runs')
#boxplot(time~n,data=serData, log="y", main="Serial wall time as a function of n",
#	   xlab="n", ylab="wall time(s)")

pdf('serial_walltime_mkl.pdf')
serData = read.csv('sermkl-1k-10k-10runs')
boxplot(time~n,data=serData, log="y", main="Serial wall time as a function of n",
	   xlab="n", ylab="wall time(s)")

#pdf('parallel_walltime1.pdf')
#parData = read.csv('par-1k-10k-10runs')
#boxplot(time.after.gather~n,data=parData, log="y", main="Parallel wall time as a function of n",
#	   xlab="n", ylab="wall time(s)")

pdf('parallel_problem_size_mkl.pdf')
parData = read.csv('parmkl-1k-10k-10runs')
boxplot(time.after.gather~n,data=parData, log="y", main="Parallel wall time as a function of n",
	   xlab="n", ylab="wall time(s)")

pdf('problem_size_speedup.pdf')
library(plyr)
serMean = ddply(serData,~n,summarize,meanTime=mean(time))
parMean = ddply(parData,~n,summarize,meanTime=mean(time.after.gather))
plot(serMean$n, serMean$meanTime/parMean$meanTime, xlab="n", ylab="speedup", main="Speedup as a function of problem size")

# changing node count stuff

#pdf('mpi_nodes_walltime.pdf')
#parNData = read.csv('par-1n-36n-2mpi-10k')
#boxplot(time.after.gather~mpi.nodes,data=parNData, log="y", main="Parallel wall time as a function of mpi processes",
#	   xlab="mpi processes", ylab="wall time(s)")

pdf('mpi_nodes_walltime_mkl.pdf')
parNData = read.csv('parmkl-1mpi-36mpi-2omp-10k-10runs')
boxplot(time.after.gather~mpi.nodes,data=parNData, log="y", main="Parallel wall time as a function of mpi processes",
	   xlab="mpi processes", ylab="wall time(s)")

pdf('mpi_nodes_speedup.pdf')
parNMean = ddply(parNData,~mpi.nodes,summarize,meanTime=mean(time.after.gather))
parNMean$mpi.nodes
rep(serMean$meanTime[10],times=6)/parNMean$meanTime
plot(parNMean$mpi.nodes, rep(serMean$meanTime[10],times=6)/parNMean$meanTime, xlab="MPI processes", ylab="speedup", main="Speedup as a function of MPI processes")
