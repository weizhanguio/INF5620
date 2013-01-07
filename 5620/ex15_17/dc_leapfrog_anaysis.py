from sympy import *
from numpy import linspace
from matplotlib.pyplot import *
a=Symbol('a')
dt=0.1
T=4
N=int(round(T/dt))
A_2=1-2*dt*a
A1=sqrt(A_2)
A2=-A1
s=A1.subs(a,1)
t=linspace(0,T,N+1)
u=linspace(0,T,N+1)
print u
c1=-1
c2=1
for n in  range(0,N+1):
	u[n]=(c1+c2*(-1)**n)*s**n	


plot(t,u)
xlabel('t')
ylabel('u')
title('C1=%g     C2=%g'   % (c1,c2))
savefig('anaysis4.jpg' )
show()






