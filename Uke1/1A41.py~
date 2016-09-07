import numpy as np

def Maxwell_Boltzmann(m,T, v):
	k=1.380649e-23
	prefac=(m/(2*np.pi*k*T))**(3/2.0)
	exponent=-0.5*(m*v**2)/(k*T)
	return prefac*np.exp(exponent)*4*np.pi*v**2


def next_step(m, T, v, delta_v):
	return (Maxwell_Boltzmann(m,T,v)+Maxwell_Boltzmann(m,T,v+delta_v))/2.0


T=15e6
n=(150e3)*1000*6.02e23
start=100
end=1000
npoints=10000
v=np.linspace(start, end, npoints)
delta_v=(end-start)/float(npoints)
integrated=0
m=1.6737e-27
for k in v:
	integrated+=next_step(m,T,k, delta_v)*delta_v
print integrated*n

