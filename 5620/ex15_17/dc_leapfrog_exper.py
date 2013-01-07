from numpy import *
from matplotlib.pyplot import *
import nose.tools as nt

def solver(I,a,b,T,dt,theta):
	dt=float(dt)
	N=int(round(T/dt))
	T = N*dt 
	u=zeros(N+1)
	t=linspace(0,T,N+1)
	u[0]=I
	u[1]=dt*(-a(t[0])*I+b(t[0]))+I
	for n in range(1,N):
		u[n+1]=(u[n-1]+2*dt*( (theta-1)*a(t[n])*u[n]+(1-theta)*b(t[n])+theta*b(t[n+1])  )) /(1+2*dt*theta*a(t[n+1]))
		#u[n+1]=u[n-1]+2*dt*(  (theta-1)*a(t[n-1])*u[n-1]-theta*a(t[n])*u[n] +(1-theta)*b(t[n-1])+theta*b(t[n]))
				
	return u,t


def test_specialcase(delta_t):       # a=1 b=1
	
	def exact_solution(t):
		return  1-exp(-t)
		
	def a(t):
		return 1.0
		
	def b(t):
		return 1.0

	theta=0;I=0;dt=delta_t
	T=4
	N=int(T/dt)
	u, t = solver(I=I, a=a, b=b, T=N*dt, dt=dt, theta=theta)
	u_e=exact_solution(t)
	return u,u_e,t
	
	
delta_t=[0.1,0.05,0.03,0.01]
for i in range(len(delta_t)):
	u,u_e,t=test_specialcase(delta_t[i])
	figure()
	plot(t,u,'r')
	plot(t,u_e)
	legend(['numerical','exact'],loc=4)
	xlabel('t')
	ylabel('u')
	title('dt=%g'   % delta_t[i])
	savefig('exper_%s.jpg' % delta_t[i])
	show()


