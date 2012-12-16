from dolfin import *
import numpy
def alpha(u):
    return  u    

rho=1.0
dt=1.0
coef=rho/phi

u = TrialFunction(V)
v = TestFunction(V)
u_k = interpolate(Constant(0.0), V)  # previous (known) u
f=Constant(0.0)
a = inner(coef*alpha(u_1)*nabla_grad(u), nabla_grad(v))*dx+u*v*dx

for n in range(1,int(T/dt)):
	
	L=v*u_1*dx+coef*f*v*dx
	u_k = interpolate(u_1, V)   #initial guess
	u = Function(V)    				 # new unknown function
	eps = 1.0           					 # error measure ||u-u_k||
	tol = 1.0E-5        					 # tolerance
	iter = 0           						 # iteration counter	
	maxiter = 50        				 # max no of iterations allowed
	while eps > tol and maxiter > iter:
		iter += 1
		solve(a == L, u)
		diff = u.vector().array() - u_k.vector().array()
		eps = numpy.linalg.norm(diff, ord=numpy.Inf)
		u_k.assign(u)                 # update for next iteration
	t=t+dt;
	u_1.assign(u)



