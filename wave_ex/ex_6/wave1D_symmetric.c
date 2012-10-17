/*
title: prove symmetry of a 1D wave problem computationally ,with u=0 at boundary 
ex6  Page:  30
*/
#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
#pragma                   message("Please input  N_x(odd),   N_t when run app")
#define   real    	      double
#define   u(i)              u[i]
#define   u_1(i)          u_1[i]
#define   u_2(i)          u_2[i]
#define   x_start        -10.0
#define   x_end        10.0 
#define   T				   1.0	                       /*time*/
#define    c                 1.0			               /* velocity*/
#define    b                 0.0	                       /*  Neumann boundary condition*/	
 

real V(real x);						                      /*initial  velocity*/
real f(real x, real t);				                      /*source term	*/
 

 
int main(int nargs, char** args)
{

	int i,n;
	int N_x, N_t;
	real dx;	 
	real dt;
	FILE *fp;


	 if (nargs>1)
       {
		   N_x = atof(args[1]);                /*number of points in x direction*/
		   N_t = atof(args[2]);                 /*number of points in time direction*/
	   }
      		printf("the number of points in N_x:%d,    N_t: %d\n",N_x, N_t);	
	
	dx=(x_end-x_start)/(N_x-1); 
	dt=T/(N_t);
	real *u=(real*)malloc((N_x+2)* sizeof(real));      
	real *u_1=(real*)malloc((N_x+2) *sizeof(real));
	real *u_2=(real*)malloc((N_x+2) *sizeof(real));
	real *temp;

	//printf("N_x=%d , N_t=%d\n",N_x, N_t);
	//printf("dx=%f, dt=%f\n",dx, dt);
	real A=pow((sqrt(c)*dt*dx),2);
			  
	if (dt>=dx/(sqrt(c)))
	     {
	   		printf("do not satisfy the requirement:c*dt/dx=1\n\n");
	     	exit(0);
	     }	
	
			
	for(i=1;i<=N_x ;i++)		                /* left boundary i=1, right boundary i=N_x*/
	    {
	      if  ( (i>=(N_x+1)/2-10)&&(i<(N_x+1)/2+10))
		       u_1(i)=1.0;       
		   else
		       u_1(i)=0.0;                           
		}	
			 				

	
	for(i=2;i<=N_x-1;i++)                                /*first time step*/	  
		{  
		u(i)=0.5*(2.0*u_1(i)+2.0*V((i-(N_x+1)/2)*dx)*dt+A*(u_1(i+1)-2.0*u_1(i)+u_1(i-1))+dt*dt*f((i-(N_x+1))/2*dx,dt));
			
		}
			
		
		temp=u_2;
		u_2=u_1;
		u_1=u;
		u=temp;
							
      
      
        for(n=1;n<=N_t;n++)
        {		 
				     for(i=2;i<=N_x-1;i++)
						{
			u(i)=2.0*u_1(i)-u_2(i)+A*(u_1(i+1)-2.0*u_1(i)+u_1(i-1))+dt*dt*f((i-(N_x+1)/2)*dx,n*dt);					 
						}		
				 
		  
			temp=u_2;
			u_2=u_1;
			u_1=u;
			u=temp;					
		}
		
		printf("store the value u at T=N_t in data.txt\n ");
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

 
 	
 real V(real x)
 	{
 	 return 1.0;	
 	}
 	
 real f(real x, real t)
 	{	 
 	 return  1.0;  
 	}
 	 
 	
