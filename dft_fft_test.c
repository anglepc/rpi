#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex.h>

double ** dft(double rex[], double imx[], int N) {
  double **X = malloc(2*sizeof(double *));
  X[0] = malloc(N*sizeof(double));
  X[1] = malloc(N*sizeof(double));

  double anglek, anglekn;
  double angle = 2.0 * M_PI / N;
  for(int k=0;k<N;k++) {
    anglek = angle * k;
    for(int n=0;n<N;n++) {
      anglekn = anglek * n;
      X[0][k] += rex[n]*cos(anglekn) + imx[n]*sin(anglekn);
      X[1][k] += -rex[n]*sin(anglekn) + imx[n]*cos(anglekn);
    }
  }

  return X;
}

complex double* dft2(complex double* x, int N) {
  complex double* X = malloc(N*sizeof(complex double*));
  malloc(N*sizeof(double));
  malloc(N*sizeof(double));

  for(int k=0;k<N;k++) {
    for(int n=0;n<N;n++) {
      X[k] += x[n] * cexp(-2.0 * M_PI * I * k * n/N);
    }
  }

  return X;
}

int* r;
double* ct;
double* st;
complex double *W;

void fft_init(char N) {
  r = malloc(N*sizeof(int));
  ct = (double *)malloc(N*sizeof(double));
  st = (double *)malloc(N*sizeof(double));
  for(int n=0;n<N;n++) {
    r[n] = 0;
    for(int j=1,k=N>>1;j<N;j<<=1,k>>=1) {
      r[n] += (n & j) ? k : 0;
    }
    ct[n] = cos(2*M_PI*n/N);
    st[n] = sin(2*M_PI*n/N);
  }

  W = malloc(N*sizeof(complex double));
  for(int k=0;k<N;k++) {
    W[k] = cexp(-2*M_PI*I*k/N);
  }
  for(int k=0;k<N;k++) {
    printf("(%+.1f %+.1fi) ", creal(W[k]), cimag(W[k]));
  }
}

complex double * fft(complex double *x, int N, int s) {
  if(N == 1) {
    return x;
  } else {
    complex double *xe = malloc(N/2*sizeof(complex double));
    complex double *xo = malloc(N/2*sizeof(comples double));
    for(int i-0;i<N;i+=2) {
      xe[i*2] = x[i];
      xo[i*2+1] = x[i+1];
    }
  }
}

complex double * fft2(complex double *x, int N int level) {
  complex double *X = malloc(N*sizeof(complex double));

// x = 0 - N-1
// y = x + :
  for(int k=0;k<N;k++) {
    for(int n=0;n<N;n++) {
      X[k] = x[r[n]] * W[k];
    }
  }
  for(int stage=1;stage<=levels;stage++) {
    for(int i=0;i<N;i++) {
    }
  }

  return X;
}

int main() {
  int N = 8; // Sample Count
  clock_t begin, end;

  double rex[] = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  double imx[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  begin = clock();
  double **X1 = dft(rex, imx, N);
  end = clock();
  printf("DFT - (%f)\n", (double)(end-begin)/CLOCKS_PER_SEC);
  for(int i=0;i<N;i++) {
    printf("  (%+.1f + %+.1fi) => (%+.1f + %+.1fi)\n", rex[i], imx[i], X1[0][i], X1[1][i]);
  }

  complex double *x = malloc(N*sizeof(complex double));
  x[0] = CMPLX(1.0, 0.0);
  x[1] = CMPLX(1.0, 0.0);
  x[2] = CMPLX(1.0, 0.0);
  x[3] = CMPLX(1.0, 0.0);
  x[4] = CMPLX(1.0, 0.0);
  x[5] = CMPLX(1.0, 0.0);
  x[6] = CMPLX(1.0, 0.0);
  x[7] = CMPLX(1.0, 0.0);

  begin = clock();
  complex double* X = dft2(x, N);
  end = clock();
  printf("DFT2 - (%f)\n", (double)(end-begin)/CLOCKS_PER_SEC);
  for(int i=0;i<N;i++) {
    printf("  (%+.1f + %+.1fi) => (%+.1f %+ .1fi) => %.1f, %.1f\n", creal(x[i]), cimag(x[i]), creal(X[i]), cimag(X[i]), cabs(X[i]), carg(X[i]));
  }

  printf("FFT Init - (%f)\n", (double)(end-begin)/CLOCKS_PER_SEC);
  fft_init(N);
  begin = clock();
  X = fft(x, N, 3);
  end = clock();
  printf("FFT - (%f)\n", (double)(end-begin)/CLOCKS_PER_SEC);
  for(int k=0;k<N;k++) {
    printf("  (%+.1f %+.1fi)\n", creal(X[k]), cimag(X[k]));
  }

  return 0;
}
