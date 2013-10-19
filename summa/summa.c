#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mkl.h"

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
   double      *p;

  // get myrow, mycol

  MPI_Comm_rank(comm_row, &mycol);
  MPI_Comm_rank(comm_col, &myrow);
  for(i=0; i<n_c[myrow]; i++)
    for(j=0; j<m_c[mycol]; j++)
      c[ldc*i+j] = beta * c[ldc*i+j];

  icurrow=0;
  icurcol=0;
  ii=0;
  jj=0;

  for(kk=0; kk<k; kk+=iwrk){
    iwrk = min(nb, m_b[icurrow]-ii);
    iwrk = min(iwrk, n_a[icurcol]-jj);

    // pack current iwrk cols of A into work1
    if(mycol == icurcol)
      memcpy(work1, a+jj, m_a[myrow]*iwrk*sizeof(double));

    // pack current iwrk rows of B into work2
    if(myrow == icurrow)
      memcpy(work2, b+(ldb*ii), iwrk*n_a[mycol]*sizeof(double));

    // broadcast work1 and work2
    RING_Bcast(work1, m_a[myrow]*iwrk, MPI_DOUBLE, icurcol, comm_row);
    RING_Bcast(work2, n_b[mycol]*iwrk, MPI_DOUBLE, icurrow, comm_col);

    // update local block
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_c[myrow], n_c[mycol],
        iwrk, alpha, work1, m_b[myrow], work2, iwrk, d_one,
        c, ldc);
    
    // update icurcol, icurrow, ii, jj
    ii += iwrk;
    jj += iwrk;
    if(jj>=n_a[icurcol]){icurcol++; jj=0;};
    if(ii>=m_b[icurrow]){icurrow++; ii=0;};
  }

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
      printf("%.20f", a[i*m+j]);
      if(j+1!=m){
        printf(", ");
      }
    }
    if(i+1!=n){
      printf(";\n");
    }
  }
  printf("];\n");
}

double *gatherMatrix(double *c, int rank, int size, int small){
  int dim = sqrt(size);
  double *tmp, *result;
  if(rank==0){
    tmp    = malloc(size*small*small*sizeof(double));
    result = malloc(size*small*small*sizeof(double));
  }
  MPI_Gather(c, small*small, MPI_DOUBLE, tmp, small*small, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(rank==0){
    // remapping of gathered result
    int ran;
    int cnt=0;
    for(ran=0; ran<size; ran++){ // for each node
      // calculate displacment
      int i_disp = small * (ran/dim); // row displacment in result
      int j_disp = small * (ran%dim); // col displacment in result
      int tmp_disp = small*dim*i_disp+j_disp; // displacment in tmp
      int i,j; // block iteration vars
      for(i=0; i<small; i++){
        for(j=0; j<small; j++){
          result[tmp_disp + i*small*dim + j] = tmp[cnt++];
        }
      }
    }
    free(tmp);
  }
  return result;
}

void scatterMatrix(double *sbuf, double *rbuf, int rank, int size, int small){
  int dim = sqrt(size);
  double *tmpSbuf;
  if(rank==0){
    tmpSbuf = malloc(size*small*small*sizeof(double));
    // remapping of gathered result
    int ran;
    int cnt=0;
    for(ran=0; ran<size; ran++){ // for each node
      // calculate displacment
      int i_disp = small * (ran/dim); // row displacment in result
      int j_disp = small * (ran%dim); // col displacment in result
      int tmp_disp = small*dim*i_disp+j_disp; // displacment in tmp
      int i,j; // block iteration vars
      for(i=0; i<small; i++)
        for(j=0; j<small; j++)
          tmpSbuf[cnt++] = sbuf[tmp_disp + i*small*dim + j];
    }
  }
  MPI_Scatter(tmpSbuf, small*small, MPI_DOUBLE, rbuf, small*small, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  free(tmpSbuf);
}

double *calcEigenVectors(int rank, int size, int small){
  double *eigs = malloc(small*small*sizeof(double));
  int dim = sqrt(size);
  int idisp = small*(rank/dim);
  int jdisp = small*(rank%dim);
  double n = dim*small;

  int i,j;
  for(i=0; i<small; i++)
    for(j=0; j<small; j++)
      eigs[i*small + j] = sqrt(2.0/(n+1.0)) * sin(((idisp+i+1.0)*(jdisp+j+1.0)*M_PI)/(n+1.0));
  return eigs;
}

double *calcTransposedEigenVectors(int rank, int size, int small){
  double *eigs = malloc(small*small*sizeof(double));
  int dim = sqrt(size);
  int idisp = small*(rank%dim);
  int jdisp = small*(rank/dim);
  double n = dim*small;

  int i,j;
  for(j=0; j<small; j++)
    for(i=0; i<small; i++)
      eigs[j*small + i] = sqrt(2.0/(n+1.0)) * sin(((jdisp+j+1.0)*(idisp+i+1.0)*M_PI)/(n+1.0));
  return eigs;
}

double *calcEigenValues(int rank, int size, int small, double h){
  double *lambda = malloc(2*small*sizeof(double));
  int dim = sqrt(size);
  int idisp = small*(rank/dim);
  int jdisp = small*(rank%dim);
  double n = dim*small;

  int j;
  for(j=0; j<small; j++){
    //lambda n j = 4/h^2 * (sin (pi*fromIntegral j/(2*(fromIntegral n+1))))^2

    //lambda[2*j] = 4.0/(h*h) * pow(sin(M_PI*(idisp+j)/(2*(n+1))),2);
    //lambda[2*j+1] = 4.0/(h*h) * pow(sin(M_PI*(jdisp+j)/(2*(n+1))),2);

    lambda[2*j] = 4.0/(h*h) * pow(sin(M_PI*(idisp+j+1)/(2*(n+1))),2);
    lambda[2*j+1] = 4.0/(h*h) * pow(sin(M_PI*(jdisp+j+1)/(2*(n+1))),2);
  }
  return lambda;
}

double *calcUmod(int size, int small, double *g, double *lambda){
  double *Umod = malloc(small*small*sizeof(double));
  int dim = sqrt(size);
  double n = dim*small;

  int i,j;
  for(i=0; i<small; i++)
    for(j=0; j<small; j++)
      Umod[i*small+j] = g[i*small+j]/(lambda[2*i]+lambda[2*j+1]);

  return Umod;
}

int main(int argc, char* argv[]){
  int rank, size;

  MPI_Init (&argc, &argv);  /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */



  int ompprocs = atoi(argv[3]);
  omp_set_num_threads(ompprocs);
  int small = atoi(argv[2]);
  int big = small*atoi(argv[1]);

  int m=big, n=big, k=big;
  double alpha=1.0;
  double beta=0.0;
  int m_a[] = {small,small,small,small,small,small,small,small,small};
  int n_a[] = {small,small,small,small,small,small,small,small,small};
  int m_b[] = {small,small,small,small,small,small,small,small,small};
  int n_b[] = {small,small,small,small,small,small,small,small,small};
  int m_c[] = {small,small,small,small,small,small,small,small,small};
  int n_c[] = {small,small,small,small,small,small,small,small,small};
  int lda = small;
  int ldb = small;
  int ldc = small;

  int nb = small;
  int dim = sqrt(size);

  double *tmp = malloc(small*small*sizeof(double));

  int i,j;

  // our space discretization step size
  double h = 1.0;

  // diagonalize the discrete laplace opperator
  double *lambda = calcEigenValues(rank, size, small, h);

  //printf("lambda%d = ", rank);
  //printMat(lambda, small, 2);
  //exit(1);

  double *eigenVectors = calcEigenVectors(rank, size, small);

  double *transEigenVectors = calcTransposedEigenVectors(rank, size, small);

  //printf("eigs%d = ", rank);
  //printMat(eigenVectors, small, small);
  //printf("\n");

  //printf("transEigs%d = ", rank);
  //printMat(transEigenVectors, small, small);
  //printf("\n");

  //exit(1);

  double *G = malloc(small*small*sizeof(double));

  //// generate initial conditions
  //srand((unsigned)time(NULL)+rank);
  //for(i=0; i<small; i++){
  //  for(j=0; j<small; j++){
  //    G[lda*i+j] = h*((double)rand()/(double)RAND_MAX);
  //  }
  //}

  // generate initial conditions
  if(rank==0){
    for(i=0; i<small; i++){
      for(j=0; j<small; j++){
        G[lda*i+j] = h*1.0;
      }
    }
  }else{
    for(i=0; i<small; i++){
      for(j=0; j<small; j++){
        G[lda*i+j] = h*0.0;
      }
    }
  }

  double *work1 = malloc(big*small*sizeof(double));
  double *work2 = malloc(big*small*sizeof(double));

  MPI_Comm comm_row, comm_col;
  MPI_Comm_split(MPI_COMM_WORLD, rank/dim, rank, &comm_row);
  MPI_Comm_split(MPI_COMM_WORLD, rank%dim, rank, &comm_col);

  double start = MPI_Wtime();
  // tmp = G * Q
  //printf("%d: tmp = G * Q\n", rank);
  pdgemm(m, n, k, nb, alpha, G, lda, eigenVectors, ldb,
    beta, tmp, ldc, m_a, n_a, m_b, n_b, m_c, n_c,
    comm_row, comm_col, work1, work2);

  double *Gmod = G;

  // Gmod = Q' * tmp
  //printf("%d: Gmod = Q' * tmp\n", rank);
  pdgemm(m, n, k, nb, alpha, transEigenVectors, lda, tmp, ldb,
    beta, Gmod, ldc, m_a, n_a, m_b, n_b, m_c, n_c,
    comm_row, comm_col, work1, work2);

  // calc Umod
  //printf("%d: calc Umod\n", rank);
  double *Umod = calcUmod(size, small, Gmod, lambda);

  // tmp = Q * Umod
  //printf("%d: tmp = Q * Umod\n", rank);
  pdgemm(m, n, k, nb, alpha, eigenVectors, lda, Umod, ldb,
    beta, tmp, ldc, m_a, n_a, m_b, n_b, m_c, n_c,
    comm_row, comm_col, work1, work2);

  double *U = Umod;
  // U = tmp * Q'
  //printf("%d: U = tmp * Q\n", rank);
  pdgemm(m, n, k, nb, alpha, tmp, lda, transEigenVectors, ldb,
    beta, U, ldc, m_a, n_a, m_b, n_b, m_c, n_c,
    comm_row, comm_col, work1, work2);

  double endMult = MPI_Wtime();

  

  double *bigU, *bigB, *bigC;
  bigU = gatherMatrix(U, rank, size, small);
  double endGather = MPI_Wtime();
  //bigB = gatherMatrix(b, rank, size, small);
  //bigC = gatherMatrix(c, rank, size, small);


  if(rank==0){
  //  //printf("A = ");
  //  //printMat(bigA, big, big);
  //  //printf("\n");
  //  //printf("B = ");
  //  //printMat(bigB, big, big);
  //  //printf("\n");
  //  //printf("C = ");
  //  //printMat(bigC, big, big);
  //  //printf("\n");
  //  FILE *fu = fopen("u.m", "w");
  //  FILE *fb = fopen("b.m", "w");
  //  FILE *fc = fopen("c.m", "w");

  //  fwrite(bigU, sizeof(double), big*big, fu);
  //  fwrite(bigB, sizeof(double), big*big, fb);
  //  fwrite(bigC, sizeof(double), big*big, fc);

  //  fclose(fu);
  //  fclose(fb);
  //  fclose(fc);
    printf("Parallel stats for %d:\n", rank);
    printf("mpi cores, omp threads, num blocks, block size, n, end mult, end gather\n");
    printf("%d, %d, %d, %d, %d, %f, %f\n", size, ompprocs, dim, small, big, endMult-start, endGather-start);
  }

  MPI_Finalize();
  return 0;
}
