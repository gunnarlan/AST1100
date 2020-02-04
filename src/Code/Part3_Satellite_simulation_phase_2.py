import numpy as np
import scipy.interpolate as inter
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem
from AST1100SolarSystemViewer import AST1100SolarSystemViewer

class Launch_simulation:

	def __init__(self, no_years, steps_per_year, seed, target_planet=6):
		self.no_years=no_years
		self.steps_per_year=steps_per_year
		self.seed=seed
		self.target_planet=target_planet
		self.mysolarsystem=AST1100SolarSystem(self.seed,hasMoons=False)
		self.major_axis=self.mysolarsystem.a
		self.mass_star=self.mysolarsystem.starMass
		self.no_planets=self.mysolarsystem.numberOfPlanets
		self.mass_of_planets=self.mysolarsystem.mass
		self.radius_planets=self.mysolarsystem.radius
		self.period_planets=self.mysolarsystem.period
		self.vx0=self.mysolarsystem.vx0
		self.vy0=self.mysolarsystem.vy0
		self.planetPos=np.load('positionsHomePlanet.npy')
		self.times=np.linspace(0, no_years, no_years*steps_per_year)
		self.positionFunction=inter.interp1d(self.times, self.planetPos)
		self.G=4*(np.pi)**2
		self.AU=1.496e11
		self.radius_planets_in_AU=1e3*self.radius_planets/float(self.AU)
		self.boost=[]
		self.errors=[]

	def compute_v_esc(self):
		radius_planet_in_AU=1e3*self.radius_planets[0]/float(self.AU)
		self.v_esc=np.sqrt(2*self.G*self.mass_of_planets[0]/float(radius_planet_in_AU))


	def compute_velocity_planets(self, t, dt):
		if t == 0:
			return np.array([self.vx0, self.vy0])
		elif t == self.no_years:
			v=(self.positionFunction(t)-self.positionFunction(t-dt))/float(dt)
			return v
		else:
			v=(self.positionFunction(t+dt)-self.positionFunction(t-dt))/float(2*dt)
			return v

	def compute_force_on_satellite(self,position, time):
		a=np.zeros(shape=(2))
		for k in xrange(self.no_planets):
			distance_to_launcher=position-self.positionFunction(time)[:,k]
			a=a-self.G*distance_to_launcher*float(self.mass_of_planets[k])/(float((np.linalg.norm(distance_to_launcher))**3))
		a=a-self.G*position*float(self.mass_star)/(float((np.linalg.norm(position))**3))
		return a

	def start_orbital_injection(self, t, position, velocity):
		planet_vel=self.compute_velocity_planets(t, 1e-8)
		relative_velocity=velocity-planet_vel[:, self.target_planet]
		position_vector=position-self.positionFunction(t)[:, self.target_planet]
		v_so=np.sqrt(self.mass_of_planets[self.target_planet]*self.G/float(np.linalg.norm(position_vector)))
		angle=np.arctan2(position_vector[1], position_vector[0])
		injection_velocity_x=v_so*np.sin(angle)-relative_velocity[0]
		injection_velocity_y=-v_so*np.cos(angle)-relative_velocity[1]
		injection_velocity=np.array([injection_velocity_x, injection_velocity_y])
		return injection_velocity

	def launch_satellite(self, N_per_year, time, t_0, initial_position, initial_v, delta_v_for_boost, delta_t_for_boost):
		ideal_dist_stable=np.linalg.norm(self.positionFunction(0)[:, self.target_planet])*np.sqrt(self.mass_of_planets[self.target_planet]/(10*float(self.mass_star)))
		ideal_dist_image=(1e3*self.radius_planets[self.target_planet]/float(self.AU))*2000.0/(float(70*np.pi/180.0))
		print "Ideal dist (stable orbit) [AU]: ", ideal_dist_stable
		print "Ideal dist (image quality) [AU]: ", ideal_dist_image
		ideal_dist=min(ideal_dist_stable, ideal_dist_image)
		print "Chosen ideal distance [AU]: ", ideal_dist

		if abs(int(N_per_year*time)-N_per_year*time) > 1e-14:
			print "Uh-oh, please ensure that the total number of timesteps is an integer!"
			sys.exit(1)

		dt=1.0/float(N_per_year)
		time_array=np.linspace(t_0, t_0+time, int(N_per_year*time))

		position_satellite=np.zeros(shape=(2, int(N_per_year*time)))
		position_satellite[:,0]=initial_position
		distance_to_target_planet=[]

		position_vector_to_planet=self.positionFunction(t_0+delta_t_for_boost)[:,self.target_planet]-initial_position
		position_vector_to_planet=position_vector_to_planet/float(np.linalg.norm(position_vector_to_planet))

		##INITIAL VELOCITY
		velocity_satellite=initial_v
		deltav=delta_v_for_boost*position_vector_to_planet
		velocity_satellite+=deltav
		already_boosted=False
		already_injected=False

		acceleration=self.compute_force_on_satellite(position_satellite[:,0], t_0)
		velocity_satellite+=acceleration*0.5*dt

		for i in xrange(0, int(N_per_year*time)-1):
			position_satellite[:,i+1]=position_satellite[:,i]+velocity_satellite*dt
			acceleration=self.compute_force_on_satellite(position_satellite[:,i+1], time_array[i+1])
			velocity_satellite=velocity_satellite+acceleration*dt

			distance_to_target_planet.append(np.linalg.norm(position_satellite[:,i+1]-self.positionFunction(time_array[i+1])[:, self.target_planet]))

			if distance_to_target_planet[-1] < ideal_dist and already_boosted==False:
				print "Yayayayayay, close enough to inject"
				injection=self.start_orbital_injection(time_array[i+1], position_satellite[:, i+1], velocity_satellite)
				print "Injection velocity [AU/year]:", np.linalg.norm(injection)
				velocity_satellite+=injection
				already_boosted=True

			if i % 1000==0:
				print "Current distance:", distance_to_target_planet[-1]

		print "Closest distance:", min(distance_to_target_planet)
		dist=np.array(distance_to_target_planet)
		min_index=np.argmin(dist)

		plt.plot(position_satellite[0,:], position_satellite[1,:])
		plt.hold('on')
		plt.plot(position_satellite[0,0], position_satellite[1,0],  marker='*')
		plt.plot(self.positionFunction(time_array)[0, self.target_planet], self.positionFunction(time_array)[1, self.target_planet])
		plt.text(position_satellite[0,0], position_satellite[1,0], "Sat")
		plt.text(self.positionFunction(t_0)[0, self.target_planet], self.positionFunction(t_0)[1, self.target_planet], "Hiffre [6]")
		plt.xlabel('Distance from star [AU]')
		plt.ylabel('Distance from star [AU]')
		plt.title('Stable orbit around target planet')
		plt.grid('on')
		plt.show()
		plt.clf()

		plt.plot(position_satellite[0,:]-self.positionFunction(time_array)[0, self.target_planet], position_satellite[1,:]-self.positionFunction(time_array)[1, self.target_planet])
		plt.hold('on')
		plt.plot(position_satellite[0,0]-self.positionFunction(t_0)[0, self.target_planet], position_satellite[1,0]-self.positionFunction(t_0)[1, self.target_planet],  marker='*')
		plt.plot(0,0,  marker='*')
		plt.text(0, 0, "Hiffre [6]")		
		plt.text(position_satellite[0,0]-self.positionFunction(t_0)[0, self.target_planet], position_satellite[1,0]-self.positionFunction(t_0)[1, self.target_planet], "Sat")
		plt.xlabel('Distance from target planet [AU]')
		plt.ylabel('Distance from target planet [AU]')
		plt.title('Stable orbit around target planet')
		plt.grid('on')




if __name__ == "__main__":
	launch=Launch_simulation(100, 20000, 6392)

	###Parameters from phase 1###
	time_start=7.09589438124
	position_start=np.array([-7.42839189752, -0.128529154053])
	velocity_start=np.array([-0.191219125118, -2.85943845417])
	delta_v=np.array(0.854) #8765
	print "Correction boost [AU/year]", delta_v

	launch.launch_satellite(N_per_year=2000000, time=0.1, t_0=time_start, initial_position=position_start, initial_v=velocity_start, delta_v_for_boost=delta_v, delta_t_for_boost=0.005435) 
