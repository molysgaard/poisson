#!/usr/bin/Rscript

library(gdata)

serData = read.csv('ser-1k-10k-20runs')
boxplot(time~n,data=serData, log="y", main="Serial wall time as a function of n",
	   xlab="n", ylab="wall time")

parData = read.csv('par-1k-10k-10runs')
boxplot(time.after.gather~n,data=parData, log="y", main="Parallel wall time as a function of n",
	   xlab="n", ylab="wall time")

library(plyr)
serMean = ddply(serData,~n,summarize,meanTime=mean(time))
parMean = ddply(parData,~n,summarize,meanTime=mean(time.after.gather))
plot(serMean$n, serMean$meanTime/parMean$meanTime, xlab="n", ylab="speedup", main="Speedup as a function of problem size")

# changing node count stuff

parNData = read.csv('par-1n-36n-2mpi-10k')
boxplot(time.after.gather~mpi.nodes,data=parNData, log="y", main="Parallel wall time as a function of mpi nodes",
	   xlab="mpi nodes", ylab="wall time")

parNMean = ddply(parNData,~mpi.nodes,summarize,meanTime=mean(time.after.gather))
parNMean$mpi.nodes
rep(serMean$meanTime[10],times=6)/parNMean$meanTime
plot(parNMean$mpi.nodes, rep(serMean$meanTime[10],times=6)/parNMean$meanTime, xlab="MPI nodes", ylab="speedup", main="Speedup as a function of MPI nodes")
