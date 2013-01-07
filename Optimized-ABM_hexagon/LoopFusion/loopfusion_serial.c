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
  real *a, *b, *y, c_j;
  real sum_a, sum_b, yp_jp1;
  clock_t start, finish;

  if (nargs>1)
    h = atof(args[1]);

  if (nargs>1)
    T = atof(args[2]);

  N = (int)(T/h+0.5);

  printf("Serial LoopFusion version---------------\n");
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

  for (j=0; j<N; j++) {
    sum_b = b[j]*f(y[0]);
    c_j = (pow(j,alpha+1)-(j-alpha)*pow(j+1,alpha))/gamma2;
    sum_a = c_j*f(y[0]);

    for (k=1; k<=j; k++) {
      sum_b += b[j-k]*f(y[k]);
      sum_a += a[j-k]*f(y[k]);
    }

    yp_jp1 = prefix + h_alpha*sum_b;
    y[j+1] = prefix + h_alpha*(sum_a + f(yp_jp1)/gamma2);
  }

  finish=clock();
  printf("y[%d]=%e",N,y[N]);
  printf("    timeusage=%g",(finish-start)/1e6);
  printf("    GFLOPs=%g\n",2.0*N*(N+1)/(finish-start)/1000.0);

  free(a);
  free(b);
  free(y);

  return 0;
}
