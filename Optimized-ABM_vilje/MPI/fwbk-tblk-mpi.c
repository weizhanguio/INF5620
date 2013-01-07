#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <mpi.h>

#define real double
#define MPI_type MPI_DOUBLE

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
  clock_t start, finish;

  int my_rank, num_procs;
  int start_k, stop_k;
  real sendbuf[4], recvbuf[4];

  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  if (nargs>1)
    h = atof(args[1]);

  if (nargs>1)
    T = atof(args[2]);

  N = (int)(T/h+0.5);

  if (my_rank==0) {
    printf("MPI Temporal-Blocking and ForwardBack combination version-----\n");
    printf("Num MPI procs=%d\n",num_procs);
    printf("T=%g, h=%g, N=%d\n", T,h,N);
  }

  a = (real*)malloc((N+1)*sizeof(real));
  b = (real*)malloc((N+1)*sizeof(real));
  y = (real*)malloc((N+1)*sizeof(real));

  h_alpha = pow(h,alpha);

  MPI_Barrier (MPI_COMM_WORLD);

  start = clock();

  for (j=N; j>=0; j--) {
    b[j] = (pow(j+1,alpha)-pow(j,alpha))/gamma1;
    a[j] = (pow(j+2,alpha+1)-2.*pow(j+1,alpha+1)+pow(j,alpha+1))/gamma2;
  }

  y[0] = y0;

  for (j=0; j<N; ) {
    if (my_rank==0) {
      sum_b = b[j]*f(y[0]);
      c_j = (pow(j,alpha+1)-(j-alpha)*pow(j+1,alpha))/gamma2;
      sum_a = c_j*f(y[0]);

      sum_b1 = b[j+1]*f(y[0]);
      c_jp1 = (pow(j+1,alpha+1)-(j+1-alpha)*pow(j+2,alpha))/gamma2;
      sum_a1 = c_jp1*f(y[0]);
    }
    else
      sum_b = sum_a = sum_b1 = sum_a1 = 0.;

    start_k = my_rank*j/num_procs+1;
    stop_k = (my_rank+1)*j/num_procs;

    for (k=start_k; k<=stop_k; k++) {
      sum_b += b[j-k]*f(y[k]);
      sum_a += a[j-k]*f(y[k]);
      sum_b1 += b[j+1-k]*f(y[k]);
      sum_a1 += a[j+1-k]*f(y[k]);
    }

    sendbuf[0] = sum_b; sendbuf[1] = sum_a;
    sendbuf[2] = sum_b1; sendbuf[3] = sum_a1;
    MPI_Allreduce (sendbuf, recvbuf, 4, MPI_type, MPI_SUM, MPI_COMM_WORLD);
    sum_b = recvbuf[0]; sum_a = recvbuf[1];
    sum_b1 = recvbuf[2]; sum_a1 = recvbuf[3];

    yp_jp1 = prefix + h_alpha*sum_b;
    y[j+1] = prefix + h_alpha*(sum_a + f(yp_jp1)/gamma2);

    yp_jp1p1 = prefix + h_alpha*(sum_b1 + b[0]*f(y[j+1]));
    y[j+2] = prefix + h_alpha*(sum_a1 + a[0]*f(y[j+1]) + f(yp_jp1p1)/gamma2);

    j += 2;

    sum_a = sum_b = sum_a1 = sum_b1 = 0.0;

    start_k = my_rank*j/num_procs+1;
    stop_k = (my_rank+1)*j/num_procs;

    for (k=stop_k; k>=start_k; k--) {
      sum_b += b[j-k]*f(y[k]);
      sum_a += a[j-k]*f(y[k]);
      sum_b1 += b[j+1-k]*f(y[k]);
      sum_a1 += a[j+1-k]*f(y[k]);
    }

    sendbuf[0] = sum_b; sendbuf[1] = sum_a;
    sendbuf[2] = sum_b1; sendbuf[3] = sum_a1;
    MPI_Allreduce (sendbuf, recvbuf, 4, MPI_type, MPI_SUM, MPI_COMM_WORLD);
    sum_b = recvbuf[0]; sum_a = recvbuf[1];
    sum_b1 = recvbuf[2]; sum_a1 = recvbuf[3];

    yp_jp1 = prefix + h_alpha*(b[j]*f(y[0]) + sum_b);
    c_j = (pow(j,alpha+1)-(j-alpha)*pow(j+1,alpha))/gamma2;
    y[j+1] = prefix + h_alpha*(c_j*f(y[0]) + sum_a + f(yp_jp1)/gamma2);

    yp_jp1p1 = prefix + h_alpha*(b[j+1]*f(y[0]) + sum_b1 + b[0]*f(y[j+1]));
    c_jp1 = (pow(j+1,alpha+1)-(j+1-alpha)*pow(j+2,alpha))/gamma2;
    y[j+2] = prefix + h_alpha*(c_jp1*f(y[0]) + sum_a1 + a[0]*f(y[j+1]) 
			       + f(yp_jp1p1)/gamma2);

    j += 2;
  }

  finish=clock();
  if (my_rank==0) {
	printf("Num MPI procs=%d",num_procs);	
    printf("	y[%d]=%e",N,y[N]);
    printf("	time usage=%g",(finish-start)/1e6);
    printf("	GFLOPs=%g\n",2.0*N*(N+1)/(finish-start)/1000.0);
  }

  free(a);
  free(b);
  free(y);

  MPI_Finalize ();

  return 0;
}
