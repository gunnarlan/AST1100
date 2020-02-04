import numpy as np
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem

N=10000

legend=['Byappo (0)', 'Domin (1)', 'Munnmon (2)', 'Pjeng (3)', 'Plaging (4)', 'Psiwe (5)', 'Hiffre (6)']
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
dt=1.0/N
for k in range(number_of_planets):
	i=1
	while x[k, i-1] < 1: 
		x[k, i]=x[k,i-1]+velocities_planets[k]*dt
		i+=1

for j in range(number_of_planets):
	noise=np.random.normal(0, 0.2, N)
	x_correct=np.logical_and(x[j,:]>=-radius_star, x[j,:]<radius_star)
	intensity[j, x_correct]=(1-((radius_planets[j])**2/(float(radius_star)**2)))
	#intensity[j,:]+=noise
	print "Minimum intensity:", min(intensity[j,:])
	print "Planet number:", j
	plt.xlabel('Position of planet [AU]')
	plt.ylabel('Intensity [%]')
	plt.title('Light curve of '+ legend[j])
	plt.grid('on')
	plt.plot(np.linspace(-1, 1, N), intensity[j,:])
	plt.show()


	
