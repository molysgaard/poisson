#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){

  int rank, count, ierr, tag;
  MPI_Comm comm;
  double sendbufA[5] = { 1, 1, 1, 1, 1 };
  double sendbufB[5] = { 2, 2, 2, 2, 2 };
  double recvbufA[5] = { 3, 3, 3, 3, 3 };
  double recvbufB[5] = { 4, 4, 4, 4, 4 };
  MPI_Status status;

  struct Matrix {
    int rows;
    int cols;
    double *data;
  };

  comm = MPI_COMM_WORLD;
  count = 5;
  tag = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &rank);

  printf("%d\n", rank);

  if(rank == 0){
    MPI_Sendrecv((void *) sendbufA, count, MPI_DOUBLE, 1, tag, (void *) recvbufA, count, MPI_DOUBLE, 1, tag, comm, &status);
  }
  else if(rank == 1){
    MPI_Sendrecv((void *) sendbufB, count, MPI_DOUBLE, 0, tag, (void *) recvbufB, count, MPI_DOUBLE, 0, tag, comm, &status);
  }
  MPI_Finalize();
}
