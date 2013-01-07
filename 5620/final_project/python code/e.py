from __future__ import division
from dolfin import *
import numpy

def b(u):
    return (u*u+1.0)
# Create mesh and define function space
mesh = Interval(40, 0,1)

V = FunctionSpace(mesh, 'Lagrange', 1)
u = TrialFunction(V)
v = TestFunction(V)

rho=1.0
t=0
dt=0.005
T=0.5
# Define variational problem for Picard iteration
u_1 = interpolate(Constant(0.0), V)  # previous (known) u
a = dt*inner(b(u_1)*nabla_grad(u), nabla_grad(v))*dx+u*v*dx

for n in range(1,int(T/dt)+1):
	t=n*dt
	f=Expression('-rho*pow(x[0],3)*0.3333 + rho*pow(x[0],2)*0.5 +  0.8889*pow(t,3)*pow(x[0],7) - 3.1111*pow(t,3)*pow(x[0],6)  +3.5*pow(t,3)*pow(x[0],5) - 1.25*pow(t,3)*pow(x[0],4) + 2*t*x[0] - t',t=t,rho=1)
	u_k = interpolate(u_1, V)   #initial guess
	L =v*u_1*dx+dt*v*f*dx
	u = Function(V)     # new unknown function
	eps = 1.0           # error measure ||u-u_k||
	tol = 1.0E-5        # tolerance
	iter = 0            # iteration counter	
	maxiter = 50        # max no of iterations allowed
	while eps > tol and maxiter > iter:
		iter += 1
		solve(a == L, u)
		diff = u.vector().array() - u_k.vector().array()
		eps = numpy.linalg.norm(diff, ord=numpy.Inf)
		#print 'iter=%d: norm=%g' % (iter, eps)
		u_k.assign(u)   # update for next iteration
	
	u_1.assign(u)



u_e=Expression('t*x[0]*x[0]*(0.5-x[0]/3)',t=T)
u_e=interpolate(u_e,V)
e = u_e.vector().array() - u.vector().array()
E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)

print "t=", T
print "dt=", dt

print "exact solution"
print u_e.vector().array()
print "     "
print "u_e-u"
print  e

#interactive()
    
