#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <gsl/gsl_cblas.h>

#define min(x,y) ( (x) < (y) ? (x) : (y) )

void printMat(double* a, int n, int m){
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

double *calcEigenVectors(int n){
  double *eigs = malloc(n*n*sizeof(double));

  int i,j;
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      eigs[i*n + j] = sqrt(2.0/(n+1.0)) * sin(((i+1.0)*(j+1.0)*M_PI)/(n+1.0));
  return eigs;
}

double *calcTransposedEigenVectors(int n){
  double *eigs = malloc(n*n*sizeof(double));

  int i,j;
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      eigs[i*n + j] = sqrt(2.0/(n+1.0)) * sin(((j+1.0)*(i+1.0)*M_PI)/(n+1.0));
  return eigs;
}

double *calcEigenValues(int n, double h){
  double *lambda = malloc(n*sizeof(double));

  int j;
  for(j=0; j<n; j++){
    lambda[j] = 4.0/(h*h) * pow(sin(M_PI*(j+1)/(2*(n+1))),2);
  }
  return lambda;
}

double *calcUmod(int n, double *g, double *lambda){
  double *Umod = malloc(n*n*sizeof(double));

  int i,j;
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      Umod[i*n+j] = g[i*n+j]/(lambda[i]+lambda[j]);

  return Umod;
}

void dmm(int n, double *A, double *B, double *C){
  int i,k,j;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      C[i*n+j] = 0;
      for(k=0; k<n; k++){
        C[i*n+j] += A[i*n+k]*B[k*n+j];
      }
    }
  }
}

int main(int argc, char* argv[]){

  int n = atoi(argv[1]);
  double alpha = 1.0;
  double beta  = 0.0;

  double *tmp = malloc(n*n*sizeof(double));

  int i,j;

  // our space discretization step size
  double h = 1.0;

  // diagonalize the discrete laplace opperator
  double *lambda = calcEigenValues(n, h);

  double *eigenVectors = calcEigenVectors(n);

  double *transEigenVectors = calcTransposedEigenVectors(n);

  double *G = malloc(n*n*sizeof(double));

  //// generate initial conditions
  //srand((unsigned)time(NULL)+rank);
  //for(i=0; i<n; i++){
  //  for(j=0; j<n; j++){
  //    G[lda*i+j] = h*((double)rand()/(double)RAND_MAX);
  //  }
  //}

  // generate initial conditions
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      G[n*i+j] = 0.0;
    }
  }
  for(i=0; i<n/6; i++){
    for(j=0; j<n/6; j++){
      G[n*i+j] = h*1.0;
    }
  }

  struct timespec before, after;
  clock_gettime(CLOCK_MONOTONIC, &before);

  // tmp = G * Q
  //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, G, n, eigenVectors, n, beta, tmp, n);
  dmm(n, G, eigenVectors, tmp);

  //double *Gmod = G;
  double *Gmod = malloc(n*n*sizeof(double));

  // Gmod = Q' * tmp
  //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, transEigenVectors, n, tmp, n, beta, Gmod, n);
  dmm(n, transEigenVectors, tmp, Gmod);

  // calc Umod
  double *Umod = calcUmod(n, Gmod, lambda);

  // tmp = Q * Umod
  //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, eigenVectors, n, Umod, n, beta, tmp, n);
  dmm(n, eigenVectors, Umod, tmp);

  double *U = Umod;
  // U = tmp * Q'
  //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, tmp, n, transEigenVectors, n, beta, U, n);
  dmm(n, tmp, eigenVectors, U);

  clock_gettime(CLOCK_MONOTONIC, &after);
  
  printf("Serial stats:\n");
  printf("n, time\n");
  printf("%d, %f\n", n, ((double) (after.tv_sec - before.tv_sec)) + (after.tv_nsec - before.tv_nsec)/10e9);

  //FILE *fu = fopen("u-seq.m", "w");
  //fwrite(U, sizeof(double), n*n, fu);
  //fclose(fu);

  return 0;
}
