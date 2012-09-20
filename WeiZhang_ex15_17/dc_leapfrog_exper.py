from numpy import *
from matplotlib.pyplot import*

def solver(I,T,dt):
	dt=float(dt)
	N=int(round(T/dt))
	u=zeros(N+1)
	u_e=zeros(N+1)
	t=linspace(0,T,N+1)
	u[0]=I
	u[1]=(1-2*dt)*u[0]+2*dt
	for n in range(0,N):
		u[n+1]=(1-2*dt)*u[n-1]+2*dt    #numerical solution
	for n in range(0,N+1):             #exact solution
		u_e[n]=1-exp(-t[n])
	
	return u, u_e,t,
	
u,u_e,t=solver(0,5,0.15)
plot(t,u_e)
plot(t,u,'r')
xlabel('t')
ylabel('u')	
legend(['exact','numerical'])
show()
