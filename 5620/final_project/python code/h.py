from dolfin import *
import numpy
import matplotlib.pyplot
def b(u):
	beta=2.0
	return (beta*u*u+1.0)


# Create mesh and define function space
mesh = UnitSquare(10, 10)
V = FunctionSpace(mesh, 'Lagrange', 2)
u = TrialFunction(V)
v = TestFunction(V)

delta=0.7**2;
u0=Expression('exp(-(x[0]*x[0]+x[1]*x[1])*0.5/delta)',delta=delta)
t=0
dt=0.02
T=0.2
# Define variational problem for Picard iteration
u_1 = interpolate(u0, V)  # previous (known) u
a = dt*inner(b(u_1)*nabla_grad(u), nabla_grad(v))*dx+u*v*dx


for n in range(1,int(T/dt)+1):
	u_k = interpolate(u_1, V)   #initial guess
	L = v*u_1*dx
	print n
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
	#t=t+dt;
	u_1.assign(u)


plot(u)


interactive()
    
