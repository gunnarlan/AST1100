import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem

N=10000

mysolarsystem=AST1100SolarSystem(6392)
radius_of_star=mysolarsystem.starRadius
number_of_planets=mysolarsystem.numberOfPlanets
initial_velocity_x=mysolarsystem.vx0
initial_velocity_y=mysolarsystem.vy0
velocities_planets=np.zeros(number_of_planets)
AU=1.5*1e11
radius_planets=1e3*(mysolarsystem.radius)/float(AU)
radius_star=1e3*(mysolarsystem.starRadius)/float(AU)
for k in range(number_of_planets):
	vel_vec=[initial_velocity_x, initial_velocity_y]
	velocities_planets[k]=np.linalg.norm(vel_vec)
	

x=np.ones(shape=[number_of_planets, N])
intensity=np.ones(shape=[number_of_planets, N])
x=-x
dt=1/float(N)
for i in range(1,N):
	x[:, i]=x[:,i-1]+velocities_planets*k*dt


for j in range(0,number_of_planets-1):
	x_correct=np.logical_and(x[j,:]>-radius_of_star, x[j,:]<radius_of_star)
	intensity[j, x_correct]=(1-((radius_planets[j])**2/(float(radius_star)**2)))
	plt.plot(x[j,:], intensity[j,:])
	plt.hold('on')

plt.show()
	
