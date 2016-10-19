##Wish to arrive at planet 6
##REMEMBER TO EXPRESS VELOCITY IN AU/YEAR!!

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
		#self.mysolarsystem.printinfo()
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
		orbital_injection=0

	def launch_satellite(self, N_per_year, time, max_close,  t_0=0, reduced_dt=1000, reduced_dt_target=100, time_to_close_planet=1.35):
		radius_planet_in_AU=1e3*self.radius_planets[0]/float(self.AU)
		#t=self.compute_start_time()
		#print "Ideal start time: ", t
		ideal_dist=np.linalg.norm(self.positionFunction(0)[:, self.target_planet])*np.sqrt(self.mass_of_planets[self.target_planet]/(10*float(self.mass_star)))

		if abs(int(N_per_year*time)-N_per_year*time) > 1e-14:
			print "Uh-oh, please ensure that the total number of timesteps is an integer!"
			sys.exit(1)

		distance_to_target_planet=[]
		number_of_years=0
		dt=1.0/float(N_per_year)
		dt_red=dt/float(reduced_dt)
		dt_red_target=dt/float(reduced_dt_target)
		time_array=np.linspace(t_0, t_0+time, int(N_per_year*time))
		self.compute_v_esc()

		position_satellite=np.zeros(shape=(2, int(N_per_year*time)))
		vel_planet=self.compute_velocity_planets(t_0, 1e-10)
		vel_home=vel_planet[:,0]
		position_satellite[:,0]=self.positionFunction(t_0)[:,0]+radius_planet_in_AU*(vel_home/float(np.linalg.norm(vel_home)))

		##INITIAL VELOCITY
		velocity_satellite=vel_home
		deltav=self.compute_needed_boost_velocity(t_0)
		print "Delta v for first boost is:", deltav-1.3134 #deltav - 0.7
		velocity_satellite+=(vel_home/(float((np.linalg.norm(vel_home)))))*(deltav-1.3134) #deltav-0.7
		already_boosted=False
		already_injected=False

		acceleration=self.compute_force_on_satellite(position_satellite[:,0], t_0)
		velocity_satellite+=0.5*acceleration*dt_red
		elapsed_time=0
		temp_pos=position_satellite[:, 0]

		for k in xrange(int(reduced_dt*N_per_year*max_close)+1):
				temp_pos+=velocity_satellite*dt_red
				elapsed_time+=dt_red
				acceleration=self.compute_force_on_satellite(temp_pos, t_0+elapsed_time)
				velocity_satellite+=acceleration*dt_red
				if k % reduced_dt == 0 and k != 0:
					position_satellite[:,int(k/reduced_dt)]=temp_pos
				distance_home=abs(np.linalg.norm(temp_pos-self.positionFunction(t_0+elapsed_time)[:,0]))
				if distance_home < radius_planet_in_AU-1e-10:
					print "Uh-oh CRASH!"
		print "Distance when initiating larger time steps:", distance_home

		velocity_satellite-=acceleration*0.5*dt_red
		velocity_satellite+=acceleration*0.5*dt
		steps_taken=int(N_per_year*max_close)

		for i in xrange(steps_taken, int(N_per_year*time)-1):
			position_satellite[:,i+1]=position_satellite[:,i]+velocity_satellite*dt
			acceleration=self.compute_force_on_satellite(position_satellite[:,i+1], time_array[i+1])
			velocity_satellite=velocity_satellite+acceleration*dt

			distance_to_target_planet.append(np.linalg.norm(position_satellite[:,i+1]-self.positionFunction(time_array[i+1])[:, self.target_planet]))

			if distance_to_target_planet[-1] < 0.01:
				print "Close enough to boost"
				iteration_when_closest=i+1
				break

			if i % 1000==0:
				print "Current distance:", distance_to_target_planet[-1]

			if i % int(N_per_year) == 0 and i != 0:
				number_of_years +=1
				print "Done with: "+str(number_of_years)+" year(s)"

		# velocity_satellite-=acceleration*0.5*dt
		# velocity_satellite+=acceleration*0.5*dt_red_target
		# elapsed_time=0
		# temp_pos=position_satellite[:, int(N_per_year*time_to_close_planet)]
		#
		# boost_direction=self.positionFunction(t_0+time_to_close_planet+0.01)[:, self.target_planet]
		# boost_magnitude=0.8
		# velocity_satellite+=((boost_direction)/np.linalg.norm(boost_direction))*boost_magnitude
		#
		# for j in xrange(int(N_per_year*time_to_close_planet), reduced_dt_target*int(N_per_year*time_remaining)-1):
		# 		temp_pos+=velocity_satellite*dt_red_target
		# 		elapsed_time+=dt_red_target
		# 		acceleration=self.compute_force_on_satellite(temp_pos, t_0+time_to_close_planet+elapsed_time)
		# 		velocity_satellite+=acceleration*dt_red
		# 		if j % reduced_dt_target == 0 and j != 0:
		# 			position_satellite[:,N_per_year*close_to_planet+int(1.0*j/reduced_dt_target)]=temp_pos
		# 		distance_target=abs(np.linalg.norm(temp_pos-self.positionFunction(t_0+time_to_close_planet+elapsed_time)[:,self.target_planet]))
		# 		if distance_target < ideal_dist:
		# 			print 'Woooop woop, we are close!'
		# 			sys.exit(0)
		#velocity_satellite-=acceleration*dt
		print "Number of crashes:", len(self.errors)
		print "Closest approach", min(distance_to_target_planet)
		print "Velocity-x:", velocity_satellite[0]
		print "Velocity-y:", velocity_satellite[1]
		print "Position-x:", position_satellite[0, iteration_when_closest]
		print "Position-y:", position_satellite[1, iteration_when_closest]
		print "Distance to target planet here:", np.linalg.norm(position_satellite[:,iteration_when_closest]-self.positionFunction(time_array[iteration_when_closest])[:, self.target_planet])
		print "Time this was achieved at:", time_array[iteration_when_closest]
		plt.plot(position_satellite[0,:], position_satellite[1,:])
		plt.hold('on')
		plt.plot(position_satellite[0,0], position_satellite[1,0],  marker='*')
		plt.text(position_satellite[0,0], position_satellite[1,0], "Sat")
		plt.hold('on')
		for k in [0,6]:
			plt.plot(self.positionFunction(self.times)[0,k], self.positionFunction(self.times)[1,k])
			plt.hold('on')
			plt.plot(self.positionFunction(t_0)[0,k], self.positionFunction(t_0)[1,k],  marker='*')
			plt.plot(self.positionFunction(t_0+time)[0,k], self.positionFunction(t_0+time)[1,k],  marker='*')
			plt.text(self.positionFunction(t_0)[0,k], self.positionFunction(t_0)[1,k], str(k))
		plt.plot(0,0, marker='*')
		plt.show()






#140 000 steps too much


launch=Launch_simulation(100, 20000, 6392)

launch.launch_satellite(200000, 4.0, max_close=1.0/100, t_0=3.5909) #1.4 #3.7451

#Chose 3.0/100 because we then are 0.1 AU away
