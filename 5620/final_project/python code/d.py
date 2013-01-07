from dolfin import *
import numpy
def a(u):
    return 1.0


# Create mesh and define function space
mesh = UnitSquare(20, 20)
V = FunctionSpace(mesh, 'Lagrange', 1)
u = TrialFunction(V)
v = TestFunction(V)

u0=Expression('cos(x[0]*pi)')
t=0
dt=0.005
T=0.5
# Define variational problem for Picard iteration
u_1 = interpolate(u0, V)  # previous (known) u
u_k = interpolate(Constant(0.0), V) 
a = dt*inner(nabla_grad(u), nabla_grad(v))*dx+u*v*dx
L = v*u_1*dx

for n in range(1,int(T/dt)):
	u_k = interpolate(u_1, V)   #initial guess
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
	t=t+dt;
	u_1.assign(u)



u_e=Expression('exp(-pi*pi*t)*cos(pi*x[0])', t=T)
u_e=interpolate(u_e,V)
e = u_e.vector().array() - u.vector().array()
E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)
print E/dt

    
