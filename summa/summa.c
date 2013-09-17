#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define A(i,j) (a[j*lda + i])
#define B(i,j) (b[j*ldb + i])
#define C(i,j) (c[j*ldc + i])

#define min(x,y) ( (x) < (y) ? (x) : (y) )

int i_one=1; //used for constant passed to blac call
double d_one=1.0,
       d_zero=0.0;

void pdgemm(m, n, k, nb, alpha, a, lda, b, ldb,
     beta, c, ldc, m_a, n_a, m_b, n_b, m_c, n_c,
     comm_row, comm_col, work1, work2)
  int m,n,k, // global matrix dims
      nb,  // panel width
      m_a[], n_a[], // dim of blocks of A
      m_b[], n_b[], // dim of blocks of A
      m_c[], n_c[], // dim of blocks of A
      lda, ldb, ldc; // leading dim of local arrays that hold local portions of matrices A, B, C
  double *a, *b, *c, // arrays that hold local parts of A, B, C
         alpha, beta, // multiplication constants
         *work1, *work2; // work arrays
  MPI_Comm comm_row, // communicator for this row of nodes
           comm_col; //  communicator for this col of nodes
{
  int myrow, mycol, // my row and col index
      nprow, npcol, // number of ned rows and cols
      i, j, kk, iwrk, // misc. index vars
      icurrow, icurcol, // index of row and col that hold current row and col, resp., for rank-1 update
      ii, jj; // local index (on icurrow and icurcol, resp.) of row and col for rank-1 update
  //double *temp, // temp pointer used in pdgemm_abt
   double      *p;

  // get myrow, mycol

  MPI_Comm_rank(comm_row, &mycol);
  MPI_Comm_rank(comm_col, &myrow);
  for(j=0; j<n_c[mycol]; j++)
    for(i=0; i<m_c[myrow]; i++)
      C(i,j) = beta * C(i,j);
  icurrow=0;
  icurcol=0;
  ii=0;
  jj=0;
  // malloc temp space for summation
  //temp = (double *) malloc(m_c[myrow]*nb*sizeof(double));

  for(kk=0; kk<k; kk+=iwrk){
    iwrk = min(nb, m_b[icurrow]-ii);
    iwrk = min(iwrk, n_a[icurcol]-jj);

    // pack current iwrk cols of A into work1
    if(mycol == icurcol)
      dlacpy_("General", &m_a[myrow], &iwrk, &A(0,jj), &lda, work1, &m_a[myrow]);

    // pack current iwrk rows of B into work2
    if(myrow == icurrow)
      dlacpy_("General", &iwrk, &n_b[mycol], &B(ii,0), &ldb, work2, &iwrk);

    // broadcast work1 and work2
    RING_Bcast(work1, m_a[myrow]*iwrk, MPI_DOUBLE, icurcol, comm_row);
    RING_Bcast(work2, n_b[mycol]*iwrk, MPI_DOUBLE, icurrow, comm_col);

    // update local block
    //dgemm_("No transpose", "No transpose", &m_c[myrow], &n_c[mycol],
    //    &iwrk, &alpha, work1, &m_b[myrow], work2, &iwrk, &d_one,
    //    c, &ldc);
    int i,j,k;
    for(i=0; i<m_c[myrow]; i++){
      for(j=0; j<n_c[mycol]; j++){
        c[n_c[mycol]*i+j] = 0;
        for(k=0; k<iwrk; k++){
          c[n_c[mycol]*i+j] += a[iwrk*i+k] * b[iwrk*k+j];
        }
      }
    }
    
    
    // update icurcol, icurrow, ii, jj
    ii += iwrk;
    jj += iwrk;
    if(jj>=n_a[icurcol]){icurcol++; jj=0;};
    if(ii>=m_b[icurrow]){icurrow++; ii=0;};
  }

  //free(temp);
}

RING_Bcast(double *buf, int count, MPI_Datatype type, int root, MPI_Comm comm){
  int me, np;
  MPI_Status status;

  MPI_Comm_rank(comm, &me);
  MPI_Comm_size(comm, &np);
  if(me != root)
    MPI_Recv(buf, count, type, (me-1+np)%np, MPI_ANY_TAG, comm, &status);
  if( (me+1)%np != root)
    MPI_Send(buf, count, type, (me+1)%np, 0, comm);
}

printMat(double* a, int n, int m){
  printf("[");
  int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<m; j++){
      printf("%f", a[i*m+j]);
      if(j+1!=m){
        printf(", ");
      }
    }
    if(i+1!=n){
      printf(";\n");
    }
  }
  printf("]\n");
}

//double *gatherMatrix(double *c, int rank, int size, int small){
//  int dim = sqrt(size);
//  double *tmp, *result;
//  if(rank==0){
//    tmp    = malloc(size*small*small*sizeof(double));
//    result = malloc(size*small*small*sizeof(double));
//  }
//  MPI_Gather(c, small*small, MPI_DOUBLE, tmp, small*small, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//  if(rank==0){
//    // remapping of gathered result
//    int ran;
//    int cnt=0;
//    for(ran=0; ran<size; ran++){ // for each node
//      // calculate displacment
//      int i_disp = small * (ran/dim); // row displacment in result
//      int j_disp = small * (ran%dim); // col displacment in result
//      int tmp_disp = small*dim*i_disp+j_disp; // displacment in tmp
//      int i,j; // block iteration vars
//      for(i=0; i<small; i++){
//        for(j=0; j<small; j++){
//          result[cnt++] = tmp[tmp_disp + i*small*dim + j];
//        }
//      }
//    }
//    free(tmp);
//  }
//  return result;
//}

int main(int argc, char* argv[]){
  int rank, size;

  MPI_Init (&argc, &argv);  /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */

  int big = 4;
  int small = 2;

  int m=big, n=big, k=big;
  double alpha=1.0;
  double beta=0.0;
  int m_a[] = {small,small,small,small};
  int n_a[] = {small,small,small,small};
  int m_b[] = {small,small,small,small};
  int n_b[] = {small,small,small,small};
  int m_c[] = {small,small,small,small};
  int n_c[] = {small,small,small,small};
  int lda = small;
  int ldb = small;
  int ldc = small;

  int nb = 2;
  int dim = sqrt(size);

  double *a = malloc(small*small*sizeof(double));
  double *b = malloc(small*small*sizeof(double));
  double *c = malloc(small*small*sizeof(double));

  int i,j;
  srand((unsigned)time(NULL)+rank);
  for(i=0; i<small; i++){
    for(j=0; j<small; j++){
      A(i,j) = ((double)rand()/(double)RAND_MAX);
      B(i,j) = ((double)rand()/(double)RAND_MAX);
    }
  }

  double *work1 = malloc(big*small*sizeof(double));
  double *work2 = malloc(big*small*sizeof(double));

  MPI_Comm comm_row, comm_col;
  MPI_Comm_split(MPI_COMM_WORLD, rank/dim, rank, &comm_row);
  MPI_Comm_split(MPI_COMM_WORLD, rank%dim, rank, &comm_col);

  double start = MPI_Wtime();
  pdgemm(m, n, k, nb, alpha, a, lda, b, ldb,
    beta, c, ldc, m_a, n_a, m_b, n_b, m_c, n_c,
    comm_row, comm_col, work1, work2);
  double endMult = MPI_Wtime();

  //double *bigA, *bigB, *bigC;
  //bigA = gatherMatrix(a, rank, size, small);
  //bigB = gatherMatrix(b, rank, size, small);
  //bigC = gatherMatrix(c, rank, size, small);

  //if(rank==0){
  //  printf("A = ");
  //  printMat(bigA, big, big);
  //  printf("\n");
  //  printf("B = ");
  //  printMat(bigB, big, big);
  //  printf("\n");
  //  printf("C = ");
  //  printMat(bigC, big, big);
  //  printf("\n");
  //}
  //double endGather = MPI_Wtime();

  //printf("Elapsed on %d: %f, %f, %f\n", rank, endMult-start, endGather-start);

  printf("A%d =", rank);
  printMat(a, small, small);
  printf("\n");

  printf("B%d =", rank);
  printMat(b, small, small);
  printf("\n");

  printf("C%d =", rank);
  printMat(c, small, small);
  printf("\n");

  MPI_Finalize();
  return 0;
}
