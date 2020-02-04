##This file contains functions taken from:
##http://stackoverflow.com/questions/16158339/plotting-a-sphere-in-python-for-an-orbital-trajectory


import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem
from Part6_Atmosphere_Model import Planet_Atmosphere
import scipy.interpolate as inter

class Orbital_Simulations():
	def __init__(self, seed, target_planet=6):
		self.seed=seed
		self.target_planet=target_planet
		self.mysolarsystem=AST1100SolarSystem(self.seed,hasMoons=False)
		self.rho_0=self.mysolarsystem.rho0[self.target_planet]
		self.M_sun_in_kg=1.989*1e30
		self.mass_of_planet_in_kg=self.M_sun_in_kg*self.mysolarsystem.mass[self.target_planet]
		self.radius_planet_in_m=1e3*self.mysolarsystem.radius[self.target_planet]
		self.G=6.674*1e-11
		self.mu=28.0101
		self.Atmosphere=Planet_Atmosphere(self.seed, self.mu, gamma=1.4)
		self.A=0.3
		self.mass_lander=90
		self.drag_coefficient=1.0
		self.drag_factor=0.5*self.drag_coefficient*self.A/float(self.mass_lander)
		self.omega=np.array([0,0,1])*(2.0*np.pi)/(24*60*60*self.mysolarsystem.period[self.target_planet]) #In m/s
		self.r_isothermal=self.Atmosphere.compute_r_T_2()

	def v_circ(self, position):
		theta=np.arctan2(position[1], position[0])
		new_v_direc=np.array([-np.sin(theta), np.cos(theta),0])
		new_v_magnitude=np.sqrt(self.G*self.mass_of_planet_in_kg/float(np.linalg.norm(position)))
		return new_v_magnitude*new_v_direc	

	def parachute_size(self, desired_v_t):
		return (2.0*self.G*self.mass_of_planet_in_kg*self.mass_lander)/(self.rho_0*self.radius_planet_in_m**2*desired_v_t**2)


	def compute_angle(self, position, time):
		r=np.linalg.norm(position)
		phi=np.arctan2(position[1], position[0])
		theta=np.arccos(position[2]/float(r))
		phi=phi-np.linalg.norm(self.omega)*time
		theta_in_deg=theta*(180.0/np.pi)
		phi_in_deg=phi*(180.0/np.pi)
		return theta_in_deg, phi_in_deg
		


	def compute_force_on_satellite(self,position, velocity):
		a=np.zeros(shape=(3))
		r=np.linalg.norm(position)
		rho=self.Atmosphere.compute_atmosphere_analytic(self.r_isothermal, r)
		v_atmosphere=np.cross(self.omega, position)
		relative_velocity=velocity-v_atmosphere
		relative_speed=np.linalg.norm(relative_velocity)
		drag_force=-self.drag_factor*rho*relative_speed*(relative_velocity)
		gravity=-self.G*position*float(self.mass_of_planet_in_kg)/(r**3)
		a=gravity+drag_force
		if np.linalg.norm(a)*self.mass_lander > 2500:
			print "Uh-oh, large force"
			print np.linalg.norm(a)*self.mass_lander
			#sys.exit(1)
		return a

	def simulate_satellite_orbit(self, N_per_second, seconds, initial_position, initial_velocity,t_0=0):
		distance_to_target_planet=np.zeros(N_per_second*seconds)
		distance_to_target_planet[0]=np.linalg.norm(initial_position)-self.radius_planet_in_m
		v_circ=self.v_circ(initial_position)
		parachute_size=self.parachute_size(2.8)
		print "Size of parachute [m]:", parachute_size

		if abs(int(N_per_second*seconds)-N_per_second*seconds) > 1e-14:
			print "Uh-oh, please ensure that the total number of timesteps is an integer!"
			sys.exit(1)
		dt=1.0/float(N_per_second)
		time_array=np.linspace(t_0, t_0+seconds, int(N_per_second*seconds))

		position_satellite=np.zeros(shape=(3, int(N_per_second*seconds)))

		already_deployed=False
		distance_before_launch=N_per_second*seconds*0.1

		fudge_factor=0.026

		##INITIAL CONDITIONS
		position_satellite[:,0]=initial_position
		velocity_satellite=initial_velocity-fudge_factor*initial_velocity  
		print "Boost needed relative to satellite:", velocity_satellite-initial_velocity

		acceleration=self.compute_force_on_satellite(position_satellite[:,0], velocity_satellite)
		velocity_satellite+=acceleration*(dt/2.0)
		self.booster_on=False
		already_circulated=False
		initial_distance=np.linalg.norm(initial_position)
		previous_step=0
		initial_theta, initial_phi=self.compute_angle(initial_position, t_0)
		print "Initial theta:", initial_theta
		print "Initial phi:", initial_phi
		
		prev_radial_vel=radial_vel=np.dot(velocity_satellite, position_satellite[:,0]/float(np.linalg.norm(position_satellite[:, 0])))
		for i in xrange(0, int(N_per_second*seconds)-1):
			position_satellite[:,i+1]=position_satellite[:,i]+velocity_satellite*dt
			v_temp=velocity_satellite+acceleration*(dt/2.0)
			acceleration=self.compute_force_on_satellite(position_satellite[:,i+1], v_temp)
			velocity_satellite=velocity_satellite+acceleration*dt
			r=np.linalg.norm(position_satellite[:, i+1])

			distance_to_target_planet[i+1]=r-self.radius_planet_in_m

			radial_vel=np.dot(velocity_satellite, position_satellite[:,i+1]/float(r))

			if radial_vel <= 0 and radial_vel_prev >=0 and already_circulated==False and i > N_per_second*seconds*0.0001:
				already_circulated=True
				print "Managed to circulate one round!"

			radial_vel_prev=radial_vel
			
			if already_deployed==False and (already_circulated==True and distance_to_target_planet[i+1] < 230*1e3):
				print "Deploying parachute!"
				print "At time [s]:", time_array[i+1]
				self.A=parachute_size
				print "Parachute size:", self.A
				self.drag_factor=0.5*self.drag_coefficient*self.A/float(self.mass_lander)
				already_deployed=True

			if distance_to_target_planet[i+1]< 0.2:
				print "Hit planet with velocity [m/s]", np.dot(velocity_satellite, position_satellite[:,i+1]/float(np.linalg.norm(position_satellite[:, i+1])))
				print "At time:", time_array[i+1]
				theta, phi=self.compute_angle(position_satellite[:, i+1], time_array[i+1])
				print "At angle theta:", theta
				print "At angle phi:", phi
				break


			if i % 100000==0:
				print "Done with:", float(i)/(N_per_second*seconds)
				print "Current force:", np.linalg.norm(acceleration)*self.mass_lander
				print "Current radial velocity:", np.dot(velocity_satellite, position_satellite[:,i+1]/float(np.linalg.norm(position_satellite[:, i+1])))
				print "Current distance:", (distance_to_target_planet[i+1])/1000.0
				print already_circulated



		###PLOTTING
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		phi = np.linspace(0, 2 * np.pi, 100)
		theta = np.linspace(0, np.pi, 100)
		xm = self.radius_planet_in_m * np.outer(np.cos(phi), np.sin(theta))
		ym = self.radius_planet_in_m * np.outer(np.sin(phi), np.sin(theta))
		zm = self.radius_planet_in_m * np.outer(np.ones(np.size(phi)), np.cos(theta))
		ax.plot_surface(xm, ym, zm)
		plt.hold('on')
		ax.plot(position_satellite[0,:], position_satellite[1,:], position_satellite[2,:], 'r--')
		ax.set_xlabel('x [m]')
		ax.set_ylabel('y [m]')
		ax.set_zlabel('z [m]')
		plt.title('Simulation of HAL9001 landing on Hiffre')
		plt.show()
		plt.plot(position_satellite[2,:],time_array/float(N_per_second*60))
		ax.set_xlabel('y [m]')
		ax.set_ylabel('Time [s]')
		plt.title('Displacement in the y-direction')
		plt.grid('on')
		plt.show()
		plt.plot(position_satellite[1,:], position_satellite[2,:])
		ax.set_xlabel('x [m]')
		ax.set_ylabel('y [m]')
		plt.grid('on')
		plt.title('xy-projection of path of HAL9001')
		plt.show()



if __name__ == "__main__":
	initial_pos_in_m=np.array([-1288663.74547901, -1946677.741875,    -819435.27243301])
	initial_vel_in_m_s=np.array([  -412.05536476,  -622.4639918,   2126.00486516])

	seconds_in_hour=60*60
	t_0=191000

	launch=Orbital_Simulations(6392)
	launch.simulate_satellite_orbit(N_per_second=200, seconds=10*seconds_in_hour, initial_position=initial_pos_in_m, initial_velocity=initial_vel_in_m_s, t_0=t_0)


