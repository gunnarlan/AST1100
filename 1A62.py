import numpy as np
import random

np.random.seed(1)

class Gas_Simulation:
	def __init__(self, L, T, N, M, t, timesteps, g, debug=False):
		self.L=L
		self.T=T
		self.N=int(N)
		self.M=M
		self.t=t
		self.timesteps=timesteps
		self.g=g

		self.k=1.38064852e-23
		self.m=2*1.6737236e-27
		self.sigma=np.sqrt(self.k*self.T/self.m)
		self.debug = debug
		self.dt=self.t/float(self.timesteps)

	def draw_particles(self):
		self.x=np.random.uniform(0, self.L, size=(self.N,3))
		self.v=np.random.normal(0, self.sigma, size=(self.N,3))
		return self.x, self.v

	def test_correctness(self):
		mean_Ek=np.mean(0.5*self.m*(np.linalg.norm(self.v, axis=1))**2)
		correct_Ek=(3/2.0)*self.k*self.T
		error_Ek = abs(mean_Ek-correct_Ek)/mean_Ek
		v_abs_mean=np.mean(np.linalg.norm(self.v, axis=1))
		correct_v_abs=np.sqrt(8*self.k*self.T/(self.m*np.pi))
		error_v=abs(v_abs_mean-correct_v_abs)/v_abs_mean
		return error_Ek, error_v

	def compute_pressure(self):
		momentum=0
		for k in range(self.timesteps):
			self.x=self.x+self.v*self.dt
			for i in range(self.x.shape[1]):
				particles_passed_right = self.x[:,i]>= self.L
				if i == 2:
					if self.debug:
						print np.sum(self.v[particles_passed_right, i])
					momentum+=2.0*self.m*(np.sum(self.v[particles_passed_right, i]))					
				self.x[particles_passed_right,i]=self.L
				self.v[particles_passed_right,i]=-self.v[particles_passed_right,i]
				particles_passed_left=self.x[:,i] <= 0
				self.x[particles_passed_left,i]=0
				self.v[particles_passed_left,i]=-self.v[particles_passed_left, i]
			if k % 100 == 0 and k!=0:
				print k
		force=momentum/(float(self.t))
		pressure=force/(float(self.L**2))
		pressure_expected=((self.N)/float(self.L**3))*self.k*self.T
		percentage_error= abs(pressure_expected-pressure)/float(pressure)
		return pressure, pressure_expected, percentage_error


	def compute_momentum_out(self):
		self.momentum_out=0
		self.number_of_particles_escaped=0
		lower_boundary_hole=self.L/4.0
		upper_boundary_hole=3*self.L/4.0
		for k in range(self.timesteps):
			self.x=self.x+self.v*self.dt
			for i in range(self.x.shape[1]):
				if i == 2:
					particles_through_hole=np.logical_and(self.x[:,2]<=0, np.logical_and(self.x[:,0]>=lower_boundary_hole, self.x[:,0]<=upper_boundary_hole), np.logical_and(self.x[:,1]>=lower_boundary_hole, self.x[:,1]<=upper_boundary_hole))
					self.momentum_out+=np.sum(self.v[particles_through_hole,2])*self.m
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
		return self.number_of_particles_escaped, self.momentum_out

	def compute_speed_gain(self, v0):
		final_velocity=v0-float(self.momentum_out)/(self.M)
		return final_velocity
			
	def compute_boxes_needed(self, v0, vf, time):
			number_of_time_steps=time/float(self.dt)
			delta_v_fuel=(self.momentum_out)/(float(number_of_particles_escaped)*self.m)
			delta_v_satelite=v0
			for i in range(number_of_time_steps):
				while delta_v_satelite <= vf:
					delta_v_satelite+=(self.number_of_particles_escaped*self.m*delta_v_fuel)/(self.M-float(self.m*i)*self.number_of_partices_escaped)
				if delta_v_satelite > vf
					break
			return number_of_boxes_needed+number_of_boxes_to_overcome_gravity
			
					
				 
		
	






simulate=Gas_Simulation(L=1e-6, T=1e4, N=1e5, M=1000, t=1e-9, timesteps=1000, g=9.81)
x_drawn, v_drawn=simulate.draw_particles()
#pressure, expected_pressure, percentage_error=simulate.compute_pressure()
#print pressure
#print expected_pressure
#print percentage_error
number_of_particles, momentum=simulate.compute_momentum_out()
velocity=simulate.compute_boxes_needed(0, 11000, 20)
print velocity

