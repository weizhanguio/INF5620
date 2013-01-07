#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define real double

#define alpha 0.5
#define gamma1 0.88622692
#define gamma2 1.329340388

#define y0 1.0
#define prefix 1.0

#define f(y) (-y)

int main(int nargs, char** args)
{
  real T=5.0, h=0.005, h_alpha;
  int j, k, N;
  real *a, *b, *y, c_j, c_jp1;
  real sum_a, sum_b, yp_jp1;
  real sum_a1, sum_b1, yp_jp1p1;
  real sum_a2, sum_b2, yp_jp1p2;
  clock_t start, finish;

  if (nargs>1)
    h = atof(args[1]);

  if (nargs>1)
    T = atof(args[2]);

  N = (int)(T/h+0.5);

  printf("Serial Temporal-Blocking version2---------------\n");
  printf("T=%g, h=%g, N=%d\n", T,h,N);

  a = (real*)malloc((N+1)*sizeof(real));
  b = (real*)malloc((N+1)*sizeof(real));
  y = (real*)malloc((N+1)*sizeof(real));

  h_alpha = pow(h,alpha);

  start = clock();

  for (j=N; j>=0; j--) {
    b[j] = (pow(j+1,alpha)-pow(j,alpha))/gamma1;
    a[j] = (pow(j+2,alpha+1)-2.*pow(j+1,alpha+1)+pow(j,alpha+1))/gamma2;
  }

  y[0] = y0;

  for (j=0; j<N; j+=2) {
    sum_b = b[j]*f(y[0]);
    c_j = (pow(j,alpha+1)-(j-alpha)*pow(j+1,alpha))/gamma2;
    sum_a = c_j*f(y[0]);

    sum_b1 = b[j+1]*f(y[0])+b[j]*f(y[1]);
    c_jp1 = (pow(j+1,alpha+1)-(j+1-alpha)*pow(j+2,alpha))/gamma2;
    sum_a1 = c_jp1*f(y[0])+a[j]*f(y[1]);

    for (k=1; k<j; k++) {
      sum_b += b[j-k]*f(y[k]);
      sum_a += a[j-k]*f(y[k]);
      sum_b1 += b[j-k]*f(y[k+1]);
      sum_a1 += a[j-k]*f(y[k+1]);
    }

    yp_jp1 = prefix + h_alpha*(sum_b + b[0]*f(y[j]));
    y[j+1] = prefix + h_alpha*(sum_a + a[0]*f(y[j]) + f(yp_jp1)/gamma2);

    yp_jp1p1 = prefix + h_alpha*(sum_b1 + b[0]*f(y[j+1]));
    y[j+2] = prefix + h_alpha*(sum_a1 + a[0]*f(y[j+1]) + f(yp_jp1p1)/gamma2);
  }

  finish=clock();
  printf("y[%d]=%e\n",N,y[N]);
  printf("time usage=%g\n",(finish-start)/1e6);
  printf("GFLOPs=%g\n",2.0*N*(N+1)/(finish-start)/1000.0);

  free(a);
  free(b);
  free(y);

  return 0;
}
