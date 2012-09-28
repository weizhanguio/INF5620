from numpy import *
from matplotlib.pyplot import*
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


def test_linear_solution():
	"""
    Test problem where u=c*t+I is the exact solution, to bereproduced by any relevant method
	"""
	def exact_solution(t):
		return c*t+I
		
	def a(t):
		return t**0.5
		
	def b(t):
		return c+a(t)*exact_solution(t)

	theta=0.4;I=0.1;dt=0.1;c=-0.5
	T=4
	N=int(T/dt)
	u, t = solver(I=I, a=a, b=b, T=N*dt, dt=dt, theta=theta)
	u_e=exact_solution(t)

	difference=abs(u_e-u).max()
	nt.assert_almost_equal(difference,0,places=14)
	
	

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
	
	sum=0.0
	for i in range(len(u)):
		sum+=(u[i]-u_e[i])**2
	E=sqrt(dt*sum)	
	return E
	
	
	
	
if __name__ == '__main__':
    test_linear_solution()
   


delta_t=[0.2,0.1,0.05,0.03,0.02,0.01]		#different dt
m=len(delta_t)
rate=zeros(m)
for i in range(m):
	rate[i]=test_specialcase(delta_t[i])
	
print "the convergence rate is  " 	
for i in range(m-1):
	converg=(log(rate[i]/rate[i+1])/log(delta_t[i]/delta_t[i+1]))
	print '%g' % converg






