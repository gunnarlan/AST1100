import numpy as np
import random
from AST1100SolarSystem import AST1100SolarSystem
from AST1100SolarSystemViewer import AST1100SolarSystemViewer

class Gas_Simulation:
	def __init__(self, L, T, N, t, timesteps,debug=False):
		np.random.seed=39232
		self.L=L
		self.T=T
		self.N=int(N)
		self.t=t
		self.timesteps=timesteps

		self.k=1.38064852e-23
		self.m_H2=2.01588 #g/mol
		self.Avogadro=6.022140e23
		self.m_molecule=self.m_H2/(float(self.Avogadro)*1000)
		self.sigma=np.sqrt(self.k*self.T/float(self.m_molecule))
		self.debug = debug
		self.dt=self.t/float(self.timesteps)
		self.x=np.random.uniform(0, self.L, size=(self.N,3))
		self.v=np.random.normal(0, self.sigma, size=(self.N,3))

	def draw_particles(self):
		self.x=np.random.uniform(0, self.L, size=(self.N,3))
		self.v=np.random.normal(0, self.sigma, size=(self.N,3))
		return self.x, self.v

	def test_correctness_of_draw(self):
		self.draw_particles()
		mean_Ek=np.mean(0.5*self.m_molecule*(np.linalg.norm(self.v, axis=1))**2)
		print "Computed kinetic energy:", mean_Ek
		correct_Ek=(3/2.0)*self.k*self.T
		print "Correct kinetic energy:", correct_Ek
		error_Ek = abs(mean_Ek-correct_Ek)/mean_Ek
		v_abs_mean=np.mean(np.linalg.norm(self.v, axis=1))
		correct_v_abs=np.sqrt(8*self.k*self.T/(self.m_molecule*np.pi))
		error_v=abs(v_abs_mean-correct_v_abs)/v_abs_mean
		return error_Ek, error_v

	def test_correctness_of_pressure(self):
		self.draw_particles()
		force=0
		percentage_error=[]
		pressure=[]
		for l in range(1):
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
					print k
			#force=momentum/(float(self.t))
			#print force
			pressure_iteration=force/(float(self.L**2)*(float(self.timesteps)))
			pressure.append(pressure_iteration)
			pressure_expected=((self.N)/float(self.L**3))*self.k*self.T
			percentage_error.append((pressure_expected-pressure_iteration)/float(pressure_expected))
			force=0
		percentage_error=np.array(percentage_error)
		pressure=np.array(pressure)
		print percentage_error
		print pressure
		return np.mean(pressure), pressure_expected, np.mean(percentage_error)

	def compute_momentum_out(self):
		self.momentum_out=0
		self.number_of_particles_escaped=0
		lower_boundary_hole=self.L/4.0
		upper_boundary_hole=3*self.L/4.0
		for k in range(self.timesteps):
			self.x=self.x+self.v*self.dt
			for i in range(self.x.shape[1]):
				if i == 2:
					particles_through_hole=np.logical_and(self.x[:,2]<=0, np.logical_and(np.logical_and(self.x[:,0]>=lower_boundary_hole, self.x[:,0]<=upper_boundary_hole), np.logical_and(self.x[:,1]>=lower_boundary_hole, self.x[:,1]<=upper_boundary_hole)))
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
				print k
		mass_per_second=(self.number_of_particles_escaped)*self.m_molecule/(float(self.t))
		force_rocket=-self.momentum_out/float(self.t)
		particles_per_second=self.number_of_particles_escaped/(float(self.t))

		print "Momentum from a single box:", self.momentum_out
		print "Number of particles escaping per second from one box:" , particles_per_second
		print "Mass per time from one box:", mass_per_second
		print "Force on rocket from single box:" ,force_rocket
		return self.momentum_out, mass_per_second, force_rocket, particles_per_second





class Launch_simulation:
	def __init__(self, momentum_out, mass_per_second, force_rocket, particles_per_second, M_sat, M_fuel,seed, L_box, no_particles, temperature,debug=False):
		self.mass_per_second=mass_per_second
		self.momentum_out=momentum_out
		self.force_rocket=force_rocket
		self.particles_per_second=particles_per_second
		self.M_sat=M_sat
		self.M_fuel=M_fuel
		self.seed=seed
		self.L_box=L_box
		self.no_particles=no_particles
		self.temperature=temperature
		self.debug=debug
		np.random.seed=39232

		self.G=6.673e-11
		self.solar_mass=1.989e30
		self.AU=1.496e11
		self.year=365*24*60*60
		self.k=1.38064852e-23

	def compute_escape_velocity_and_g(self):
		mySolarSystem=AST1100SolarSystem(self.seed)
		print mySolarSystem.mass[0]
		print mySolarSystem.radius[0]
		mass_of_my_planet=self.solar_mass*mySolarSystem.mass[0]
		radius_of_my_planet=1e3*mySolarSystem.radius[0]
		print "Mass of planet:", mass_of_my_planet
		print "Radius of planet:", radius_of_my_planet
		self.v_esc = np.sqrt(2*self.G*mass_of_my_planet/float(radius_of_my_planet))
		self.g=self.G*(mass_of_my_planet)/(float(radius_of_my_planet)**2)
		print "Escape velocity:", self.v_esc
		print "g:" , self.g
		return self.v_esc, self.g


	def compute_boxes_needed_for_escape_ignore_fuel(self, time,N):
		time=time*60
		steps=time/float(N)
		#Analytic
		self.Number_of_boxes_ignore_fuel=np.ceil(self.M_sat*self.v_esc/(self.force_rocket*float(time)))
		print "Number of boxes needed if ignoring fuel (analytic):", self.Number_of_boxes_ignore_fuel

		#Numeric
		v=0
		for k in range(N):
			a=(self.force_rocket/float(self.M_sat))
			v+=a*steps
		self.Number_of_boxes_ignore_fuel_numeric=(self.v_esc)/float(v)
		print "Number of boxes needed if ignoring fuel (numeric):", self.Number_of_boxes_ignore_fuel_numeric


	def compute_boxes_needed_for_escape_with_fuel(self, time, N):
		time=time*60
		steps=float(time)/N
		v=0
		#Compute velocity reached by a single box
		for k in range(N):
			mass=self.M_sat+self.M_fuel-k*steps*(self.mass_per_second)
			a=(self.force_rocket/float(mass))
			v+=a*steps
		self.Number_of_boxes_with_fuel=np.ceil(float(self.v_esc)/v)
		print "Number of boxes, not ignoring fuel:", self.Number_of_boxes_with_fuel

	def compute_boost_ignore_mass_loss(self,step_size, delta_v, v0):
		self.delta_v_ignore_mass_loss=delta_v
		F=self.Number_of_boxes_ignore_fuel*self.force_rocket
		fuel_per_second=self.Number_of_boxes_ignore_fuel*self.mass_per_second
		print "Fuel per second:", fuel_per_second
		mass=self.M_sat
		v=v0
		steps_taken=0
		self.fuel_used_ignore_mass_loss=0
		self.time_taken_ignore_mass_loss=0
		while np.linalg.norm(v) < np.linalg.norm(v0+delta_v):
			self.time_taken_ignore_mass_loss+=step_size
			v+=step_size*F/float(mass)
			self.fuel_used_ignore_mass_loss+=fuel_per_second*step_size
		print "Fuel used (ignoring mass loss):", self.fuel_used_ignore_mass_loss

	def test_consistency(self):
		system=AST1100SolarSystem(self.seed)
		system.massNeededCheck(float(self.Number_of_boxes_ignore_fuel), float(self.delta_v_ignore_mass_loss), float(self.force_rocket), float(self.particles_per_second), float(self.fuel_used_ignore_mass_loss))


	def compute_boost_with_mass_loss(self, step_size, deltav, v0, m_fuel): #Ignoring gravity
		self.step_size_with_mass_loss=step_size
		self.deltav=deltav
		F=self.Number_of_boxes_ignore_fuel_numeric*self.force_rocket
		print "Force from boxes:", F
		self.v_with_mass_loss=[v0]
		time_taken=0
		i=0
		while np.linalg.norm(self.v_with_mass_loss[i]) < np.linalg.norm(deltav):
			i+=1
			time_taken+=step_size
			mass=self.M_sat+m_fuel-i*step_size*(self.mass_per_second)*self.Number_of_boxes_ignore_fuel_numeric
			self.v_with_mass_loss.append(self.v_with_mass_loss[i-1]+(float(F)/(mass))*step_size)
		self.mass_lost=time_taken*(self.mass_per_second)*self.Number_of_boxes_ignore_fuel_numeric
		if self.debug:
			print self.v_with_mass_loss

		return time_taken, self.mass_lost

	def compute_boost_analytic(self, v_0, force_per_box, mass_out_per_box, delta_v):
		exponent=(delta_v*mass_out_per_box)/float(force_per_box)
		return 1100*(np.exp(exponent)-1)

	def create_movie(self):
		system=AST1100SolarSystemViewer(self.seed)
		system.escapeVelMovie(self.v_with_mass_loss, self.step_size_with_mass_loss)
		print "Done!"




L_box=1e-6
N=1e5
T=1e4

simulate=Gas_Simulation(L=L_box, T=T, N=N, t=1e-9, timesteps=1000)
simulate.draw_particles()
simulate.test_correctness_of_draw()
momentum_out, mass_per_second, force_from_rocket, particles_per_second=simulate.compute_momentum_out()
pressure, pressure_expected, percentage_error=simulate.test_correctness_of_pressure()
print "Analytical pressure:", pressure_expected
print "Calculated pressure:", pressure
print "Percentage error:", percentage_error

launch=Launch_simulation(momentum_out, mass_per_second, force_from_rocket, particles_per_second, M_sat=1100, M_fuel=17199.5129366, seed=6392, L_box=L_box, no_particles=N,temperature=T)
v_esc, g=launch.compute_escape_velocity_and_g()
delta_v=4.753691274*(1.496e11)/float(365*24*60*60)
print "Delta v:", delta_v
delta_m=launch.compute_boost_analytic(v_0=0, force_per_box=force_from_rocket, mass_out_per_box=mass_per_second, delta_v=delta_v+0.1*delta_v)
print "Mass needed (analytic):" , delta_m
launch.compute_boxes_needed_for_escape_ignore_fuel(20, 10000)
#launch.compute_boxes_needed_for_escape_with_fuel(20, 100000)
#launch.compute_boost_ignore_mass_loss(step_size=1e-4, delta_v=1000, v0=0)
time, mass=launch.compute_boost_with_mass_loss(step_size=0.001, deltav=delta_v+0.1*delta_v, v0=0, m_fuel=23000)
print mass
#launch.test_consistency()
#launch.create_movie()
