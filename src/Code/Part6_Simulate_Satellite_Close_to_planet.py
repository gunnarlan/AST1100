import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem
from Part6_Atmosphere_Model import Planet_Atmosphere
import scipy.interpolate as inter

class Orbital_Simulations():
	def __init__(self, seed, target_planet=6, atmosphere_tolerance=1e-13):
		self.seed=seed
		self.target_planet=target_planet
		self.mysolarsystem=AST1100SolarSystem(self.seed,hasMoons=False)
		self.M_sun_in_kg=1.989*1e30
		self.mass_of_planet_in_kg=self.M_sun_in_kg*self.mysolarsystem.mass[self.target_planet]
		self.radius_planet_in_m=1e3*self.mysolarsystem.radius[self.target_planet]
		self.G=6.674*1e-11
		self.mu=28.0101
		self.Atmosphere=Planet_Atmosphere(self.seed, self.mu, gamma=1.4)
		self.A=6 
		self.drag_coefficient=1.0
		self.mass_satellite=1100.0
		self.drag_factor=0.5*self.drag_coefficient*self.A/float(self.mass_satellite)
		self.omega=np.array([0,0,1])*(2.0*np.pi)/(24*60*60*self.mysolarsystem.period[self.target_planet])
		self.r_isothermal=self.Atmosphere.compute_r_T_2()
		self.r_critical=self.Atmosphere.compute_critical_r(atmosphere_tolerance)
		print "Radius:", self.radius_planet_in_m/1000.0
		print "Critical r (km)", (self.r_critical-self.radius_planet_in_m)/1000.0

	def compute_first_boost_hohmann(self, position):
		r2=np.linalg.norm(position)
		r1=self.r_critical
		mu=self.G*self.mass_of_planet_in_kg
		deltav=np.sqrt(mu/float(r2))*(1-np.sqrt(2.0*r1/(r1+r2)))
		transfer_time=np.pi*np.sqrt((r1+r2)**3/(8.0*mu))
		print "Deltav:", deltav
		return deltav, transfer_time


	def compute_second_boost_hohmann(self, position):
		r2=np.linalg.norm(position)
		r1=self.r_critical
		mu=self.G*self.mass_of_planet_in_kg
		deltav=np.sqrt(mu/float(r1))*(-1+np.sqrt(2.0*r2/(r1+r2)))
		transfer_time=np.pi*np.sqrt((r1+r2)**3/(8.0*mu))
		print "Deltav:", deltav
		print "Transfer time:", transfer_time
		return deltav,transfer_time


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
		ratio_forces=np.linalg.norm(gravity)/float(np.linalg.norm(drag_force))
		if ratio_forces < 1000:
			print "You are too close, the atmosphere has a big effect here"
			print (np.linalg.norm(position)-self.radius_planet_in_m)/1000.0
			sys.exit(1)
		if np.linalg.norm(a) > 250000:
			print "Uh-oh, your lander died"
			sys.exit(1)
		return a

	def v_circ(self, position, velocity):
		theta=np.arctan2(position[1], position[0])
		new_v_direc=np.array([-np.sin(theta), np.cos(theta),0])
		new_v_magnitude=np.sqrt(self.G*self.mass_of_planet_in_kg/float(np.linalg.norm(position)))
		return new_v_magnitude*new_v_direc-velocity

	def v_circ_polar(self, position, velocity):
		new_v_direc=np.array([0.0,0.0,1.0])
		new_v_magnitude=np.sqrt(self.G*self.mass_of_planet_in_kg/float(np.linalg.norm(position)))
		return new_v_magnitude*new_v_direc-velocity		

	def simulate_satellite_orbit(self, N_per_second, seconds, initial_position, initial_velocity,t_0=0):
		already_boosted=False

		distance_to_target_planet=[]
		deltav_1, transfer_time=self.compute_first_boost_hohmann(initial_position)
		deltav_2=self.compute_second_boost_hohmann(initial_position)
		v_circ=self.v_circ(initial_position, initial_velocity)

		if abs(int(N_per_second*seconds)-N_per_second*seconds) > 1e-14:
			print "Uh-oh, please ensure that the total number of timesteps is an integer!"
			sys.exit(1)
		dt=1.0/float(N_per_second)
		time_array=np.linspace(t_0, t_0+seconds, int(N_per_second*seconds))

		position_satellite=np.zeros(shape=(3, int(N_per_second*seconds)))

		already_boosted=False

		##INITIAL CONDITIONS
		position_satellite[:,0]=initial_position
		velocity_satellite=initial_velocity
		velocity_satellite+=v_circ
		print "First boost:", v_circ
		velocity_satellite-=deltav_1*(velocity_satellite)/float(np.linalg.norm(velocity_satellite))
		print "Second boost:", -deltav_1*(velocity_satellite)/float(np.linalg.norm(velocity_satellite))


		acceleration=self.compute_force_on_satellite(position_satellite[:,0], velocity_satellite)
		velocity_satellite+=acceleration*(dt/2.0)

		for i in xrange(0, int(N_per_second*seconds)-1):
			position_satellite[:,i+1]=position_satellite[:,i]+velocity_satellite*dt
			v_temp=velocity_satellite+acceleration*(dt/2.0)
			acceleration=self.compute_force_on_satellite(position_satellite[:,i+1], v_temp)
			velocity_satellite=velocity_satellite+acceleration*dt

			distance_to_target_planet.append(np.linalg.norm(position_satellite[:,i+1]))


			if i % 100000==0:
				print "Current distance:", (distance_to_target_planet[-1]-self.radius_planet_in_m)/1000.0
				print "Done with:", float(i)/(N_per_second*seconds)
			if i*dt>transfer_time and already_boosted==False:
				deltav=self.v_circ_polar(position_satellite[:, i+1], velocity_satellite) 
				velocity_satellite += deltav
				already_boosted=True
				print "Boosting!, at ",  time_array[i]
				print "Boosting with", deltav

		print "Minimum distance to planet [km]:", (min(distance_to_target_planet)-self.radius_planet_in_m)/1000.0
		print "Maximum distance to planet [km]:", (max(distance_to_target_planet)-self.radius_planet_in_m)/1000.0
		print "Minimum distance to planet [m]:", (min(distance_to_target_planet))
		print "Maximum distance to planet [m]:", (max(distance_to_target_planet))
		print "Final position [m]:", position_satellite[:, -1]
		print "Final velocity [m/s]:'", velocity_satellite

		###PLOTTING
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.set_zlim(-6*1e7, 6*1e7)
		u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
		x=self.radius_planet_in_m*np.sin(v)*np.cos(u)
		y=self.radius_planet_in_m*np.sin(u)*np.sin(v)
		z=self.radius_planet_in_m*np.cos(v)
		u, v = np.mgrid[0:2*np.pi:10j, 0:np.pi:5j]
		x_atm=self.r_critical*np.sin(v)*np.cos(u)
		y_atm=self.r_critical*np.sin(u)*np.sin(v)
		z_atm=self.r_critical*np.cos(v)
		ax.plot_surface(x, y, z, color="b")
		plt.hold('on')
		ax.plot_wireframe(x_atm, y_atm, z_atm, color="r")
		ax.plot(position_satellite[0,:], position_satellite[1,:], position_satellite[2,:])
		ax.set_xlabel('x [m]')
		ax.set_ylabel('y [m]')
		ax.set_zlabel('z [m]')
		plt.title('Satellite orbiting planet')
		plt.show()



if __name__ == "__main__":

	[planetPos, t]=np.load('planetPositions.npy')
	positionFunction=inter.interp1d(t, planetPos)
	target_planet=6
	time_0_in_years=9.0
	dt=1e-12
	years_in_seconds=365.24*24*60*60
	AU_in_m=1.495978707*1e11

	initial_pos_in_m=np.array([34454171.6279 ,  52030864.2132 ,  0 ]) #From init
	initial_vel_in_m_s=np.array([ 344.694950787 ,  -266.905761513 ,  0])

	seconds_in_year=365*30*60*60 #Approximate value to avoid decimals in time steps.

	launch=Orbital_Simulations(6392)
	launch.simulate_satellite_orbit(N_per_second=10, seconds=0.01*seconds_in_year, initial_position=initial_pos_in_m, initial_velocity=initial_vel_in_m_s, t_0=0)
