#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <omp.h>

#define real double

#define alpha 0.5
#define gamma1 0.88622692
#define gamma2 1.329340388

#define y0 1.0
#define prefix 1.0

#define f(y) (-y)

int main(int nargs, char** args)
{
  real T=5.0, h=0.000025, h_alpha;
  int j, k, N;
  real *a, *b, *y, c_j, c_jp1;
  real sum_a, sum_b, yp_jp1;
  real sum_a1, sum_b1, yp_jp1p1;
  double start, finish;
  //omp_set_num_threads(4);
  if (nargs>1)
    h = atof(args[1]);

  if (nargs>1)
    T = atof(args[2]);

  N = (int)(T/h+0.5);

  printf("OpenMP Temporal-Blocking version---------------\n");
  printf("T=%g, h=%g, N=%d\n", T,h,N);

  a = (real*)malloc((N+1)*sizeof(real));
  b = (real*)malloc((N+1)*sizeof(real));
  y = (real*)malloc((N+1)*sizeof(real));

  h_alpha = pow(h,alpha);

#pragma omp parallel
  {
#pragma omp master
    printf("Number of threads = %d", omp_get_num_threads());
  }

  start = omp_get_wtime();

#pragma omp parallel for
  for (j=N; j>=0; j--) {
    b[j] = (pow(j+1,alpha)-pow(j,alpha))/gamma1;
    a[j] = (pow(j+2,alpha+1)-2.*pow(j+1,alpha+1)+pow(j,alpha+1))/gamma2;
  }

  y[0] = y0;
  sum_a = sum_b = sum_a1 = sum_b1 = 0.0;

#pragma omp parallel default(shared) private(j,k)
    {
  for (j=0; j<N; j+=2) {

#pragma omp for reduction(+: sum_b,sum_a,sum_b1,sum_a1)
    for (k=1; k<=j; k++) {
      sum_b += b[j-k]*f(y[k]);
      sum_a += a[j-k]*f(y[k]);
      sum_b1 += b[j+1-k]*f(y[k]);
      sum_a1 += a[j+1-k]*f(y[k]);
    }

#pragma omp single
    {
    
    yp_jp1 = prefix + h_alpha*(b[j]*f(y[0]) + sum_b);
    c_j = (pow(j,alpha+1)-(j-alpha)*pow(j+1,alpha))/gamma2;
    y[j+1] = prefix + h_alpha*(c_j*f(y[0]) + sum_a + f(yp_jp1)/gamma2);

    yp_jp1p1 = prefix + h_alpha*(b[j+1]*f(y[0]) + sum_b1 + b[0]*f(y[j+1]));
    c_jp1 = (pow(j+1,alpha+1)-(j+1-alpha)*pow(j+2,alpha))/gamma2;
    y[j+2] = prefix + h_alpha*(c_jp1*f(y[0]) + sum_a1 + a[0]*f(y[j+1])
			       + f(yp_jp1p1)/gamma2);

    sum_a = sum_b = sum_a1 = sum_b1 = 0.0;
    
    }
  }
    }  /* end of the parallel region */

  finish=omp_get_wtime();
  
  printf("	y[%d]=%e",N,y[N]);
  printf("	time usage=%g",(finish-start));
  printf("	GFLOPs=%g\n",2.0*N*(N+1)/(finish-start)/1e9);

  free(a);
  free(b);
  free(y);

  return 0;
}
