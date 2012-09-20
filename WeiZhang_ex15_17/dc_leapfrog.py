from numpy import *
from matplotlib.pyplot import*
def solver(I,T,dt):
	dt=float(dt)
	N=int(round(T/dt))
	u=zeros(N+1)
	u_e=zeros(N+1)
	t=linspace(0,T,N+1)
	sum=0.0
	u[0]=I
	u[1]=(1-2*dt)*u[0]+2*dt
	for n in range(0,N):
		u[n+1]=(1-2*dt)*u[n-1]+2*dt    #numerical solution
	for n in range(0,N+1):             #exact solution
		u_e[n]=1-exp(-t[n])
	for n in range(0,N+1):
		sum+=(u_e[n]-u[n])**2
	E=sqrt(dt*sum)
	return E





M=50
a=linspace(0.01,0.5,M)
e=zeros(len(a))
for i in range(0,M):
	e[i]=solver(0,5,a[i])
r=log(e[1]/e[0])/log(a[1]/a[0])
print r

