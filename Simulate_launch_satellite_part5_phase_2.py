##Wish to arrive at planet 6
##REMEMBER TO EXPRESS VELOCITY IN AU/YEAR!!

import numpy as np
import scipy.interpolate as inter
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem


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
		[self.planetPos, t]=np.load('planetPositions.npy')
		self.positionFunction=inter.interp1d(t, self.planetPos)
		self.G=4*(np.pi)**2
		self.AU=1.496e11
		self.radius_planets_in_AU=1e3*self.radius_planets/float(self.AU)
		self.boost=[]
		self.errors=[]

	def compute_v_esc(self):
		radius_planet_in_AU=1e3*self.radius_planets[0]/float(self.AU)
		self.v_esc=np.sqrt(2*self.G*self.mass_of_planets[0]/float(radius_planet_in_AU))
		#print "Escape velocity:", self.v_esc


	def compute_needed_boost_velocity(self, t0):
		v_plan=self.compute_velocity_planets(t0, 1e-10)
		v_home=np.linalg.norm(v_plan[:,0])
		r_target=np.linalg.norm(self.positionFunction(t0)[:,self.target_planet])
		r_home=np.linalg.norm(self.positionFunction(t0)[:,0])
		delta_v_1=v_home*(np.sqrt(2.0*r_target/(r_target+float(r_home))-1))
		mu=self.G*self.mass_of_planets[0]
		delta_v_rel_earth=np.sqrt(delta_v_1**2+2*mu/self.radius_planets_in_AU[0])
		return delta_v_rel_earth

	def compute_needed_angle(self):
		mu=self.G*self.mass_star
		a_p=self.major_axis[0]
		home_P=2*np.pi*np.sqrt(a_p**3/float(mu))
		target_P=2*np.pi*np.sqrt(self.major_axis[self.target_planet]**3/float(mu))
		target_omega=2*np.pi/float(target_P)
		angle_offset=np.pi-target_omega*home_P
		angle_offset_wiki=np.pi*(1-1.0/(2*np.sqrt(2))*np.sqrt((self.major_axis[0]/float(self.major_axis[self.target_planet])+1)**3))
		print "The correct angle offset is", angle_offset
		print "The correct angle offset from wiki is:", angle_offset_wiki
		return angle_offset_wiki

	def compute_start_time(self):
		angle_wanted=self.compute_needed_angle()
		eps=1e-3
		years=10
		number_of_steps_per_year=10000
		dt=1.0/number_of_steps_per_year
		angle_home=np.arctan2(self.positionFunction(0)[1,0], self.positionFunction(0)[0,0])
		angle_goal=np.arctan2(self.positionFunction(0)[1, self.target_planet], self.positionFunction(0)[0, self.target_planet])
		angle_diff =angle_goal-angle_home
		time_steps_passed=0
		while abs(angle_wanted-angle_diff) > eps and time_steps_passed < number_of_steps_per_year*years:
			time_steps_passed+=1
			angle_home=np.arctan2(self.positionFunction(time_steps_passed*dt)[1,0], self.positionFunction(time_steps_passed*dt)[0,0])
			angle_goal=np.arctan2(self.positionFunction(time_steps_passed*dt)[1, self.target_planet], self.positionFunction(time_steps_passed*dt)[0, self.target_planet])
			angle_diff =angle_goal-angle_home

		t_final=time_steps_passed*dt

		if abs(angle_wanted-angle_diff) <= eps:
			return t_final
		else:
			print "Failed to converge"
			sys.exit(1)

	def compute_velocity_planets(self, t, dt):
		if t == 0:
			return np.array([self.vx0, self.vy0])
		elif t == self.no_years:
			v=(self.positionFunction(t)-self.positionFunction(t-dt))/float(dt)
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
		print "Injection velocity:", injection_velocity
		return injection_velocity


	def launch_satellite(self, N_per_year, time, t_0, initial_position, initial_v, delta_v_for_boost, delta_t_for_boost):
		ideal_dist=np.linalg.norm(self.positionFunction(0)[:, self.target_planet])*np.sqrt(self.mass_of_planets[self.target_planet]/(10*float(self.mass_star)))
		print "Ideal dist is:", ideal_dist

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

		initial_dist=np.linalg.norm(position_satellite[:,0]-self.positionFunction(time_array[0])[:, self.target_planet])
		print "Initial distance:", initial_dist

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
				print "Time:", time_array[i+1]/float(N_per_year*time)+t_0
				injection=self.start_orbital_injection(time_array[i+1], position_satellite[:, i+1], velocity_satellite)
				print "Injection velocity:", np.linalg.norm(injection)
				velocity_satellite+=injection
				already_boosted=True

			if i % 1000==0:
				print "Current distance:", distance_to_target_planet[-1]

		print "Closest distance:", min(distance_to_target_planet)
		dist=np.array(distance_to_target_planet)
		min_index=np.argmin(dist)

		#plt.plot(position_satellite[0,:]-self.positionFunction(time_array)[0, self.target_planet], position_satellite[1,:]-self.positionFunction(time_array)[1, self.target_planet])
		#plt.hold('on')
		#plt.plot(position_satellite[0,0]-self.positionFunction(t_0)[0, self.target_planet], position_satellite[1,0]-self.positionFunction(t_0)[1, self.target_planet],  marker='*')
		#plt.text(position_satellite[0,0], position_satellite[1,0], "Sat")
		#plt.hold('on')

		plt.plot(position_satellite[0,:], position_satellite[1,:])
		plt.hold('on')
		plt.plot(position_satellite[0,0], position_satellite[1,0],  marker='*')
		plt.plot(position_satellite[0,min_index], position_satellite[1,min_index],  marker='*')
		plt.text(position_satellite[0,min_index], position_satellite[1,min_index], str(min_index))

		plt.plot(self.positionFunction(time_array)[0,self.target_planet], self.positionFunction(time_array)[1,self.target_planet])
		#plt.text(self.positionFunction(time_array[0])[0,self.target_planet], self.positionFunction(time_array[0])[1,self.target_planet])
		plt.show()

		#plt.plot(position_satellite[0,min_index+10], position_satellite[1,min_index+10],  marker='*')
		#plt.text(position_satellite[0,min_index+10], position_satellite[1,min_index+10], str(min_index+10))

		#plt.plot(self.positionFunction(time_array[min_index+10])[0,self.target_planet], self.positionFunction(time_array[min_index+10])[1,self.target_planet],  marker='*')
		#plt.text(self.positionFunction(time_array[min_index+10])[0,self.target_planet], self.positionFunction(time_array[min_index+10])[1,self.target_planet], str(min_index+10))




		#for k in [0,6]:
		#	plt.plot(self.positionFunction(self.times)[0,k], self.positionFunction(self.times)[1,k])
		#	plt.hold('on')
		#	plt.plot(self.positionFunction(t_0)[0,k], self.positionFunction(t_0)[1,k],  marker='*')
		#	plt.plot(self.positionFunction(t_0+time)[0,k], self.positionFunction(t_0+time)[1,k],  marker='*')
		#	plt.text(self.positionFunction(t_0)[0,k], self.positionFunction(t_0)[1,k], str(k))
		#plt.plot(0,0, marker='*')
		#plt.show()





launch=Launch_simulation(100, 20000, 6392)

#time_start=5.00039948321
#position_start=np.array([-3.62295006174, 6.31873695092])
#velocity_start=np.array([-5.10803644909, 1.82897151167])
#delta_v=np.array(4.7291)

#time_start=4.99867447705
#position_start=np.array([-3.61413746601, 6.31557959543])
#velocity_start=np.array([-5.10952263349, 1.83174958298])
#delta_v=np.array(1.4393) #1.4393

time_start=7.09589438124
position_start=np.array([-7.4144299907, -0.149121699967])
velocity_start=np.array([-0.177788924778, -2.86452352387])
delta_v=0#np.array(0.854) #8765


launch.launch_satellite(N_per_year=2000000, time=0.1, t_0=time_start, initial_position=position_start, initial_v=velocity_start, delta_v_for_boost=delta_v, delta_t_for_boost=0.05435) #0.005435

#Chose 3.0/100 because we then are 0.1 AU away
