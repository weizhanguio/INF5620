/*
title: use ghost cells to implement Neumann boundary conditions, 1D wave with 
Neumann boundary on both side
ex4  Page:  30
*/
#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
#pragma                   message("Please input  N_x,   N_t when run app")
#define   real    	      double
#define   x(i)              x[i]
#define   u(i)              u[i]
#define   u_1(i)          u_1[i]
#define   u_2(i)          u_2[i]
#define   x_start        0.0
#define   x_end        10.0 
#define   T				   1.0					  /*time*/
#define    c                 1.0                     /* velocity*/			
#define    b                 0.0	                  /*  Neumann boundary condition*/		
 
real I(real x );        					              /*initial  value of u*/
real V(real x );						                  /*initial  velocity*/
real f(real x, real t);				                  /*source term	*/
 

 
int main(int nargs, char** args)
{

	int i,n;
	int N_x, N_t;
	real dx;	 
	real dt;
	FILE *fp;


	 if (nargs>1)
       {
		   N_x = atof(args[1]);      /*number of points in x direction*/
		   N_t = atof(args[2]);        /*number of points in time direction*/
	   }
      		
	
	dx=(x_end-x_start)/(N_x-1); 
	dt=T/(N_t);
	real *u=(real*)malloc((N_x+2)* sizeof(real));      
	real *u_1=(real*)malloc((N_x+2) *sizeof(real));
	real *u_2=(real*)malloc((N_x+2) *sizeof(real));
	real *temp;

	printf("%d points in x direction, %d points in time direction\n",N_x, N_t);
	
	real A=pow((sqrt(c)*dt*dx),2);    /* square of Courant Number  */
			  
	if (dt>=dx/(sqrt(c)))
	     {
	   		printf("do not satisfy the requirement:c*dt/dx<=dx\n\n");
	     	exit(0);
	     }	
	
			
	for(i=1;i<=N_x ;i++)		                                   /* left boundary i=1, right boundary i=N_x*/
		  u_1(i)=I(i*dx);                                              /*implement initial condition*/
			
			 				
		u_1(0)=u_1(2)+2.0*b*dx;								/*implement left ghost cells */
		u_1(N_x+1)=u_1(N_x-1)+2.0*b*dx;           
	
			
	for(i=1;i<=N_x;i++)                                   /* first time step*/	  
		{  
		u(i)=0.5*(2.0*u_1(i)+2.0*V(i*dx)*dt+A*(u_1(i+1)-2.0*u_1(i)+u_1(i-1))+dt*dt*f(i*dx,dt));
			
		}
			
 			 /*update ghost points  for n=1*/
			u(0)=u(2)+2.0*b*dx;								 
			u(N_x+1)=u(N_x)+2.0*b*dx;                 
		 
		
		temp=u_2;
		u_2=u_1;
		u_1=u;
		u=temp;
							
        for(n=2;n<=N_t;n++)
        {		 
				     for(i=1;i<=N_x;i++)
						{
							u(i)=2.0*u_1(i)-u_2(i)+A*(u_1(i+1)-2.0*u_1(i)+u_1(i-1))+dt*dt*f(i*dx,n*dt);					 
						}		
				 
		  /*update ghost points */
			u(0)=u(2)+2.0*b*dx;								 
			u(N_x+1)=u(N_x)+2.0*b*dx; 
			
			      
			temp=u_2;
			u_2=u_1;
			u_1=u;
			u=temp;					
		}
		 
		 
	printf("store the value u at T=N_t in data.txt \n");	 
	fp=fopen("data.txt","w")	;
      for (i=1;i<=N_x;i++)
      	    {
      	    fprintf(fp,"%f\n",u_1(i));
            }
        
    fclose(fp);	
	
		free(u);
		free(u_1);
		free(u_2);	 
		temp=NULL;
	

}



 
 real I(real x)
 	{
 	 return sin(x);
 	}
 	
 real V(real x)
 	{
 	 return 1.0;	
 	}
 	
 real f(real x, real t)
 	{	 
 	 return  1.0;  
 	}
 	 
 	
