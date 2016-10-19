###SWITCHED TO EULER-CROMER, NEED BOTH VEL AND A
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem
from Atmosphere_Model import Planet_Atmosphere
import scipy.interpolate as inter

class Orbital_Simulations():
	def __init__(self, seed, target_planet=6, atmosphere_tolerance=1e-5):
		self.seed=seed
		self.target_planet=target_planet
		self.mysolarsystem=AST1100SolarSystem(self.seed,hasMoons=False)
		self.M_sun_in_kg=1.989*1e30
		self.mass_of_planet_in_kg=self.M_sun_in_kg*self.mysolarsystem.mass[self.target_planet]
		self.radius_planet_in_m=1e3*self.mysolarsystem.radius[self.target_planet]
		self.G=6.674*1e-11
		self.mu=28.0101
		self.Atmosphere=Planet_Atmosphere(self.seed, self.mu, gamma=1.4)
		self.A=6 #Rember to change this value
		self.drag_coefficient=1.0
		self.drag_factor=0.5*self.drag_coefficient*self.A
		self.omega=np.array([0,0,1])*(2.0*np.pi)/(24*60*60*self.mysolarsystem.period[self.target_planet])
		self.r_isothermal=self.Atmosphere.compute_r_T_2()
		self.r_critical=self.Atmosphere.find_root_closest_r()
		self.ratio_forces=1e60

	def compute_first_boost_hohmann(self, position):
		mu=self.G*self.mass_of_planet_in_kg
		r=np.linalg.norm(position)
		deltav=np.sqrt(mu/float(self.r_critical))*(1-np.sqrt(2.0*self.r_critical/(self.r_critical+r)))
		return deltav

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
		if ratio_forces < self.ratio_forces:
			self.ratio_forces=ratio_forces
		if ratio_forces < 1000:
			print "You are too close, the atmosphere has a big effect here"
			sys.exit(1)

		if np.linalg.norm(a) > 250000:
			print "Uh-oh, your lander died"
			sys.exit(1)
		return a

	def simulate_satellite_orbit(self, N_per_second, seconds, initial_position, initial_velocity,t_0=0):
		already_boosted=False
		
		distance_to_target_planet=[]
		deltav_1=self.compute_first_boost_hohmann(initial_position)
		
		if abs(int(N_per_second*seconds)-N_per_second*seconds) > 1e-14:
			print "Uh-oh, please ensure that the total number of timesteps is an integer!"
			sys.exit(1)
		dt=1.0/float(N_per_second)
		time_array=np.linspace(t_0, t_0+seconds, int(N_per_second*seconds))

		position_satellite=np.zeros(shape=(3, int(N_per_second*seconds)))


		##INITIAL CONDITIONS
		position_satellite[:,0]=initial_position
		velocity_satellite=initial_velocity

		for i in xrange(0, int(N_per_second*seconds)-1):
			acceleration=self.compute_force_on_satellite(position_satellite[:,i], velocity_satellite)
			velocity_satellite=velocity_satellite+acceleration*dt
			position_satellite[:,i+1]=position_satellite[:,i]+velocity_satellite*dt

			distance_to_target_planet.append(np.linalg.norm(position_satellite[:,i+1]))

			if i % 1000==0:
				print "Current distance:", distance_to_target_planet[-1]
				print "Done with:", float(i)/(N_per_second*seconds)
				#velocity_satellite = velocity_satellite-5*(position_satellite[:,i]/(float(np.linalg.norm(position_satellite[:,i]))))
			if distance_to_target_planet[-1] < 2367000 and already_boosted==False:
				velocity_satellite = velocity_satellite-0.1623*velocity_satellite
				already_boosted=True
				print "Boosting!"
		
		print min(distance_to_target_planet)
		print self.r_critical
        
		###PLOTTING
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		#ax.set_xlim(-self.radius_planet_in_m, self.radius_planet_in_m)
		#ax.set_ylim(-self.radius_planet_in_m, self.radius_planet_in_m)
		#ax.set_zlim(-6*1e7, 6*1e7)
		#plt.gca().set_aspect('equal', adjustable='box')
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
		#plt.savefig('Close_orbit_satellite.pdf')
		plt.show()
		print self.ratio_forces



[planetPos, t]=np.load('planetPositions.npy')
positionFunction=inter.interp1d(t, planetPos)
target_planet=6
time_0_in_years=9.0
dt=1e-12
years_in_seconds=365.24*24*60*60
AU_in_m=1.495978707*1e11

v_plan=(positionFunction(time_0_in_years+dt)-positionFunction(time_0_in_years-dt))/float(2*dt)
initial_pos_in_au=np.array([-4.42891695786 , -6.14234439191, 0])-np.array([positionFunction(time_0_in_years)[0, target_planet], positionFunction(time_0_in_years)[1, target_planet],0])
initial_vel_in_au_year=np.array([2.80237588654 , -2.14120396904,0])-np.array([v_plan[0, target_planet],v_plan[1, target_planet], 0])
#initial_vel_in_au_year=np.array([0,0,0])

time_0_in_seconds=years_in_seconds*time_0_in_years
initial_pos_in_m=np.array([ 652366.35859064, -5084567.5092116,         0.        ])
initial_vel_in_m_s=np.array([-1197.64507154,  -627.28957862,     0.        ])

seconds_in_year=300*30*60*60
#2364757.81961

launch=Orbital_Simulations(6392)
launch.simulate_satellite_orbit(N_per_second=1, seconds=0.002*seconds_in_year, initial_position=initial_pos_in_m, initial_velocity=initial_vel_in_m_s, t_0=time_0_in_seconds)
