import math
import random
import numpy as np

def Gauss(x,sigma, mu):
	prefac=1/(sigma*math.sqrt(2*math.pi))
	exponent=-((x-mu)**2)/(2.0*sigma**2)
	return prefac*math.exp(exponent)


def next_step(x, delta_x, sigma, mu):
	return (Gauss(x,sigma,mu)+Gauss(x+delta_x,sigma,mu))/2.0

sigma=2
mu=29
start=19
end=39
npoints=10000
x=np.linspace(start, end, npoints)
delta_x=(end-start)/float(npoints)
integrated=0
for k in x:
	integrated+=next_step(k, delta_x, sigma, mu)*delta_x
print integrated
