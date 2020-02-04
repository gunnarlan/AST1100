import numpy as np
import sys
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
from AST1100SolarSystemViewer import AST1100SolarSystemViewer

class Gas_Simulation:
	def __init__(self, L, T, N, t, timesteps,debug=False):
		self.L=L #Length of box
		self.T=T #Temperature of particles
		self.N=int(N) #Number of particles
		self.t=t #Simulation time
		self.timesteps=timesteps
		self.debug = debug


		self.k=1.38064852e-23 #Boltzmann
		self.m_H2=2.01588 #g/mol
		self.Avogadro=6.022142*1e23
		self.m_molecule=self.m_H2/(float(self.Avogadro)*1000)
		self.sigma=np.sqrt(self.k*self.T/float(self.m_molecule)) #Ideal gas
		self.dt=self.t/float(self.timesteps)
		self.x=np.random.uniform(0, self.L, size=(self.N,3))
		self.v=np.random.normal(0, self.sigma, size=(self.N,3))
		if self.debug==True:
			np.random.seed=39399

	def draw_particles(self):
		self.x=np.random.uniform(0, self.L, size=(self.N,3))
		self.v=np.random.normal(0, self.sigma, size=(self.N,3))

	def test_correctness_of_draw(self):
		"""
		Test function to check if the distribution behaves as expected
		"""

		print "Testing the draw"
		mean_Ek=np.mean(0.5*self.m_molecule*(np.linalg.norm(self.v, axis=1))**2)
		print "Computed kinetic energy [J]:", mean_Ek
		correct_Ek=(3/2.0)*self.k*self.T
		print "Correct kinetic energy [J]:", correct_Ek
		error_Ek = abs((mean_Ek-correct_Ek)/float(mean_Ek))
		print "Percentage error in kinetic energy [%]:", error_Ek
		v_abs_mean=np.mean(np.linalg.norm(self.v, axis=1))
		correct_v_abs=np.sqrt(8*self.k*self.T/(self.m_molecule*np.pi))
		error_v=100*abs((v_abs_mean-correct_v_abs)/float(v_abs_mean))
		print "Percentage error in average speed [%]:", error_v

	def test_correctness_of_pressure(self, number_of_iterations):
		"""
		Test function to check correct implementation of pressure calcuations
		"""

		print "Testing pressure implementation with %d iterations" %(number_of_iterations)
		force=0
		pressure=np.zeros(number_of_iterations)
		pressure_expected=((self.N)/float(self.L**3))*self.k*self.T
		for l in range(1,number_of_iterations+1):
			self.draw_particles()
			print "Iteration number:", l
			for k in range(self.timesteps):
				self.x=self.x+self.v*self.dt
				for i in range(self.x.shape[1]):
					particles_passed_right = self.x[:,i]>= self.L
					if i == 2:
						if self.debug:
							print np.sum(self.v[particles_passed_right, i])
						force+=2.0*self.m_molecule*(np.sum(self.v[particles_passed_right, i]))/float(self.dt)
					self.x[particles_passed_right,i]=self.L
					self.v[particles_passed_right,i]=-self.v[particles_passed_right,i]
					particles_passed_left=self.x[:,i] <= 0
					self.x[particles_passed_left,i]=0
					self.v[particles_passed_left,i]=-self.v[particles_passed_left, i]
				if k % 100 == 0 and k!=0:
					print "Progres of pressure calculations: " + str(100.0*k/self.timesteps)+'%'
			pressure_iteration=force/(float(self.L**2)*(float(self.timesteps)))
			pressure[l-1]=pressure_iteration
			pressure_expected=((self.N)/float(self.L**3))*self.k*self.T
			force=0
		if self.debug:
			print pressure
		mean_pressure=np.mean(pressure)
		percentage_error=100*abs((mean_pressure-pressure_expected)/float(pressure_expected))
		print "Mean pressure with %2d iterations is %.5f [Pa]" %(number_of_iterations, mean_pressure)
		print "Analytic pressure should be: "+ str(pressure_expected)+" [Pa]"
		print "Percentage error [%]:", np.mean(percentage_error)

	def compute_momentum_out(self):
		self.momentum_out=0
		self.number_of_particles_escaped=0
		lower_boundary_hole=self.L/4.0
		upper_boundary_hole=3*self.L/4.0
		for k in range(self.timesteps):
			self.x=self.x+self.v*self.dt
			for i in range(self.x.shape[1]):
				if i == 2: #Z-direction (arbtirary)
					particles_through_hole=np.logical_and(self.x[:,2]<=0, np.logical_and(np.logical_and(self.x[:,0]>=lower_boundary_hole, self.x[:,0]<=upper_boundary_hole), 							np.logical_and(self.x[:,1]>=lower_boundary_hole, self.x[:,1]<=upper_boundary_hole)))
					self.momentum_out+=np.sum(self.v[particles_through_hole,2])*self.m_molecule
					if self.debug:
						print np.sum(particles_through_hole)
					self.number_of_particles_escaped+=np.sum(particles_through_hole)
				particles_passed_right = self.x[:,i]>= self.L
				self.x[particles_passed_right,i]=self.L
				self.v[particles_passed_right,i]=-self.v[particles_passed_right,i]
				particles_passed_left=self.x[:,i] <= 0
				self.x[particles_passed_left,i]=0
				self.v[particles_passed_left,i]=-self.v[particles_passed_left, i]
			if k % 100 == 0 and k!=0:
				print "Progres of momentum calculations: " + str(100.0*k/self.timesteps)+'%'
		mass_per_second=(self.number_of_particles_escaped)*self.m_molecule/(float(self.t))
		force_rocket=-self.momentum_out/float(self.t)
		particles_per_second=self.number_of_particles_escaped/(float(self.t))
		return self.momentum_out, mass_per_second, force_rocket, particles_per_second





class Launch_simulation:
	def __init__(self, momentum_out, mass_per_second, force_rocket, particles_per_second, M_sat,seed, L_box, no_particles, temperature,debug=False):
		self.mass_per_second=mass_per_second
		self.momentum_out=momentum_out
		self.force_rocket=force_rocket
		self.particles_per_second=particles_per_second
		self.M_sat=M_sat
		self.L_box=L_box
		self.no_particles=no_particles
		self.temperature=temperature
		self.debug=debug
		self.seed=seed
		if self.debug:
			np.random.seed=39232

		self.G=6.673e-11
		self.solar_mass=1.989e30
		self.AU=1.496e11
		self.year=365.24*24*60*60
		self.k=1.38064852e-23
		self.m_H2=2.01588 #g/mol
		self.Avogadro=6.022142*1e23
		self.m_molecule=self.m_H2/(float(self.Avogadro)*1000)
		self.Number_of_boxes_ignore_fuel=0
		self.Number_of_boxes_with_fuel=0

	def compute_escape_velocity_and_g(self):
		mySolarSystem=AST1100SolarSystem(seed)
		mass_of_my_planet=self.solar_mass*mySolarSystem.mass[0]
		radius_of_my_planet=1e3*mySolarSystem.radius[0]
		self.v_esc = np.sqrt(2*self.G*mass_of_my_planet/float(radius_of_my_planet))
		self.g=self.G*(mass_of_my_planet)/(float(radius_of_my_planet)**2)
		return self.v_esc, self.g


	def compute_boxes_needed_for_escape_ignore_fuel(self, time,N):
		"""
		Here we ignore fuel loss
		"""
		v,g=self.compute_escape_velocity_and_g()
		time=time*60 #Time in minutes
		steps=time/float(N)
		#Analytic
		self.Number_of_boxes_ignore_fuel=np.ceil((self.M_sat)*self.v_esc/(self.force_rocket*float(time)))
		print "Number of boxes needed if ignoring fuel (analytic):", self.Number_of_boxes_ignore_fuel

		#Numeric
		v=0
		for k in range(N):
			a=(self.force_rocket/float(self.M_sat))
			v+=a*steps
		self.Number_of_boxes_ignore_fuel_numeric=(self.v_esc)/float(v)
		print "Number of boxes needed if ignoring fuel (numeric):", self.Number_of_boxes_ignore_fuel_numeric


	def compute_boxes_needed_for_escape_with_fuel(self, M_fuel, time, N):
		"""
		Here we account for fuel loss
		"""

		time=time*60
		steps=float(time)/int(N)
		v=0
		#Compute velocity reached by a single box
		for k in range(int(N)):
			mass=self.M_sat+M_fuel-k*steps*(self.mass_per_second)
			a=(self.force_rocket/float(mass))
			v+=a*steps
		self.Number_of_boxes_with_fuel=np.ceil(float(self.v_esc)/v)
		print "Number of boxes, not ignoring fuel:", self.Number_of_boxes_with_fuel

	def compute_boost_ignore_mass_loss(self,N, delta_v, v0=0):
		step_size=1.0/N
		if self.Number_of_boxes_ignore_fuel == 0:
			print "Error, please compute the number of boxes first"
			sys.exit(1)
		self.delta_v_ignore_mass_loss=delta_v
		F=self.Number_of_boxes_ignore_fuel*self.force_rocket
		fuel_per_second=self.Number_of_boxes_ignore_fuel*self.particles_per_second*self.m_molecule
		print "Fuel per second:", fuel_per_second
		mass=self.M_sat
		v=v0
		steps_taken=0
		self.fuel_used_ignore_mass_loss=0.0
		self.time_taken_ignore_mass_loss=0.0
		while np.linalg.norm(v) < np.linalg.norm(v0+delta_v):
			self.time_taken_ignore_mass_loss+=step_size
			v+=step_size*F/float(mass)
			self.fuel_used_ignore_mass_loss+=fuel_per_second*step_size
		print "Fuel used (ignoring mass of fuel) [kg]:", self.fuel_used_ignore_mass_loss

	def compute_boost_with_mass_loss(self, step_size, deltav, v0, m_fuel, plot=True, AU_or_m_s='AU'): #Ignoring gravity
		if self.Number_of_boxes_with_fuel == 0:
			print "Error, please compute the number of boxes first"
			sys.exit(1)
		self.step_size_with_mass_loss=step_size
		self.deltav=deltav
		F=self.Number_of_boxes_with_fuel*self.force_rocket
		self.v_with_mass_loss=[v0]
		mass=[self.M_sat+m_fuel]
		time_taken=0
		i=0
		while np.linalg.norm(self.v_with_mass_loss[i]) < np.linalg.norm(deltav):
			i+=1
			time_taken+=step_size
			mass.append(self.M_sat+m_fuel-i*step_size*(self.mass_per_second)*self.Number_of_boxes_with_fuel)
			self.v_with_mass_loss.append(self.v_with_mass_loss[i-1]+(float(F)/(mass[-1]))*step_size)
		self.mass_lost=time_taken*(self.mass_per_second)*self.Number_of_boxes_with_fuel
		if self.debug:
			print self.v_with_mass_loss

		if plot:
			if AU_or_m_s=='AU':
				self.v_with_mass_loss=self.year*np.array(self.v_with_mass_loss)/(float(self.AU))
			time_array=np.linspace(0, time_taken, len(self.v_with_mass_loss))
			plt.plot(time_array, self.v_with_mass_loss)
			plt.xlabel('Time, t [s]')
			if AU_or_m_s=='AU':
				plt.plot(time_array, np.ones(len(time_array))*self.deltav*self.year/(float(self.AU)), 'r--')
				plt.ylabel('Velocity, v [AU/year]')
			else:
				plt.plot(time_array, np.ones(len(time_array))*self.deltav*self.year(float(self.AU)), 'r--')
				plt.ylabel('Velocity, v [m/s]')
			plt.title('Velocity of rocket escaping from planet, ignoring gravity')
			plt.grid('on')
			plt.legend(['Velocity achieved by rocket', 'Required velocity'], prop={'size':10})
			plt.savefig('Rocket_vel.pdf')
			#plt.show()
			plt.clf()
			plt.plot(time_array, mass)
			plt.plot(time_array, np.ones(len(time_array))*self.M_sat, 'r--')
			plt.xlabel('Time, t [s]')
			plt.ylabel('Mass remaining, m [kg]')
			plt.title('Mass loss as a function of time for the rocket ignoring gravity')
			plt.grid('on')
			plt.legend(['Total mass remaining', 'Mass of VoyagerX'], prop={'size':10})
			plt.savefig('Rocket_mass.pdf')
			#plt.show()

		print "Mass remaining at the end of boost [kg]:", mass[-1]

		return time_taken, self.mass_lost

	def test_consistency(self):
		"""
		Test function to test reasonableness
		"""

		system=AST1100SolarSystem(self.seed)
		self.compute_boxes_needed_for_escape_ignore_fuel(20, 100000)
		self.compute_boost_ignore_mass_loss(1000, self.v_esc, 0)
		system.massNeededCheck(float(self.Number_of_boxes_ignore_fuel), self.v_esc, float(self.force_rocket), float(self.particles_per_second), float(self.fuel_used_ignore_mass_loss))

	def compute_boost_analytic(self, v_0, force_per_box, mass_out_per_box, delta_v):
		exponent=(delta_v*mass_out_per_box)/float(force_per_box)
		return self.M_sat*(np.exp(exponent)-1)

	def create_movie(self):
		system=AST1100SolarSystemViewer(self.seed)
		system.escapeVelMovie(self.v_with_mass_loss, self.step_size_with_mass_loss)
		print "Done!"




if __name__ == "__main__":
	L_box=1e-6
	No_particles=1e5
	N=1000
	T=1e4
	check_implementation=False
	seed=6392
	delta_v_AU_year=4.753691274
	delta_v_m_s=delta_v_AU_year*(1.496e11)/float(365.24*24*60*60) #From simulation of orbit
	time_to_reach_esc=20 #minutes
	safety_percentage=0.1 #10% safety margin in case of errors



	simulate=Gas_Simulation(L=L_box, T=T, N=No_particles, t=1e-9, timesteps=1000)
	if check_implementation:
		simulate.test_correctness_of_draw()
		simulate.test_correctness_of_pressure(2)
	momentum_out, mass_per_second, force_from_rocket, particles_per_second=simulate.compute_momentum_out()

	print "Force from a single box [kg m/s]:", force_from_rocket
	print "Number of particles escaping per time from one box [1/s]:" , particles_per_second
	print "Mass per time from one box [kg/s]:", mass_per_second
	print "Force on rocket from one box [J]:", force_from_rocket
	launch=Launch_simulation(momentum_out, mass_per_second, force_from_rocket, particles_per_second, M_sat=1100, seed=6392, L_box=L_box, no_particles=No_particles,temperature=T)
	print "----------------------------------------"
	print "        Starting the rocket             "
	print "----------------------------------------"
	boost=delta_v_m_s+safety_percentage*delta_v_m_s
	print "Trying to achieve boost of [AU/year]: ", delta_v_AU_year+safety_percentage*delta_v_AU_year
	v_esc, g=launch.compute_escape_velocity_and_g()
	delta_m=launch.compute_boost_analytic(v_0=0, force_per_box=force_from_rocket, mass_out_per_box=mass_per_second, delta_v=boost) #10% insurance
	print "Mass needed (analytic):" , delta_m
	launch.compute_boxes_needed_for_escape_with_fuel(delta_m, time_to_reach_esc, N)
	time, mass=launch.compute_boost_with_mass_loss(step_size=0.001, deltav=boost, v0=0, m_fuel=delta_m)
	print "Mass needed (numeric) [kg]: ", mass
	print "Relative difference in mass calculation [%]", 100*abs((delta_m-mass)/float(delta_m))

	if check_implementation:
		launch.test_consistency()
		launch.create_movie()

