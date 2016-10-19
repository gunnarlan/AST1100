import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem
from AST1100SolarSystemViewer import AST1100SolarSystemViewer


class Calculate_planet_orbits:

	def __init__(self, seed):
		self.seed=seed
		self.mysolarsystem=AST1100SolarSystem(self.seed,hasMoons=False)
		self.mass_star=self.mysolarsystem.starMass
		self.no_planets=self.mysolarsystem.numberOfPlanets
		self.mass_of_planets=self.mysolarsystem.mass
		self.eccentricity=self.mysolarsystem.e
		self.major_axis=self.mysolarsystem.a
		self.radius_of_star_in_AU=self.mysolarsystem.starRadius*1.0/(1.5*1e11)
		self.radius_of_planets_in_AU=self.mysolarsystem.radius*1.0/(1.5*1e11)
		self.x0=self.mysolarsystem.x0
		self.y0=self.mysolarsystem.y0
		self.vx0=self.mysolarsystem.vx0
		self.vy0=self.mysolarsystem.vy0
		self.mysolarsystem.printinfo()
		self.G=4*(np.pi)**2

	def compute_force_one_body(self,position, m_object_producing_force):
		acceleration=np.zeros_like(position)
		if len(position.shape) < 1:
			print "Error, the acceleration function requires at least a 2D array"
			sys.exit(1)
		elif len(position.shape)==1: #Only one planet
			distance=np.linalg.norm(position)
			a=-self.G*m_object_producing_force/float(distance**3)
			a_vec=a*position
			acceleration=a_vec

		elif len(position.shape) == 2: #Multiple planets
			for i in xrange(len(position[0])):
				distance=np.linalg.norm(position[:,i])
				a=-self.G*m_object_producing_force/float(distance**3)
				a_vec=a*position[:,i]
				acceleration[:,i]=a_vec
		else:
			print "Not understood datatype"
			sys.exit(1)
		return acceleration

	def compute_force_two_bodies(self, position, m_planets):
		if len(position.shape) <2:
			print "Error, to compute a two-body problem you need at least two bodies!"
			sys.exit(1)
		elif len(position.shape) ==2:
			acceleration=np.zeros_like(position)
			for i in xrange(len(m_planets)):
				a=np.zeros(2)
				if i != len(m_planets)-1:
					distance=position[:,i]-position[:,-1]
					a=-(self.G*m_planets[-1]/(float(np.linalg.norm(distance)**3)))*distance
				else:
					for j in xrange(len(m_planets)-1):
						distance=position[:,j]-position[:,-1]
						a=a+(self.G*m_planets[j]/(float(np.linalg.norm(distance)**3)))*distance
				acceleration[:, i]=a
		else:
			print "Not understood datatype"
			sys.exit(1)
		return acceleration


	def compute_force_many_bodies(self, position, m_planets):
		if len(position.shape) <2:
			print "Error, to compute a multi-body problem you need at least two bodies!"
			sys.exit(1)
		elif len(position.shape) ==2:
			acceleration=np.zeros_like(position)
			for i in xrange(len(m_planets)):
				a=np.zeros(2)
				for k in xrange(len(m_planets)):
					if k != i:
						distance=position[k]-position[i]
						a+=(self.G*m_planets[k]/(float(np.linalg.norm(distance))**3))*distance
				acceleration[i]=a
		else:
			print "Not understood datatype"
			sys.exit(1)
		return acceleration



	def  compute_orbits_one_body(self, time, N, debug=False, tol=0.01, upper=1, lower=1e-15):


		#NB Time in years
		dt=time/float(N)
		position_of_planets=np.zeros(shape=(2, self.no_planets, N))
		position_of_planets[0,:,0]=self.x0
		position_of_planets[1,:,0]=self.y0

		v=np.zeros(shape=(2,self.no_planets))
		v[0,:]=self.vx0
		v[1,:]=self.vy0

		if debug:
			error=1
			pos_ref=[]
			npoints=np.log10(upper)-np.log10(lower)
			h=[(10)**-j for j in xrange(int(np.ceil(npoints))+1)]
			required_h=np.zeros(self.no_planets)
			k=0
			for l in xrange(self.no_planets):
				acceleration=self.compute_force_one_body(position_of_planets[:,l,0], self.mass_star)
				while error > tol and k < npoints:
					pos_ref.append(position_of_planets[:,l,0]+v[1,:]*h[k]+0.5*acceleration*h[k]**2)
					if k != 0:
						norm_current=np.linalg.norm(pos_ref[k])
						norm_previous=np.linalg.norm(pos_ref[k-1])
						error=abs((norm_current-norm_previous)/float(norm_previous))
					k+=1
				if k>= npoints:
					print "Not succesful, try smaller h"
					sys.exit(1)
				else:
					required_h[l]=(h[k])
			print "Succesful, smallest h needed is:", required_h.min()




		#Implement leapfrog method
		i=0
		no_years=0
		acceleration=self.compute_force_one_body(position_of_planets[:,:,0], self.mass_star)
		v+=0.5*acceleration*dt
		for i in xrange(1,N):
			position_of_planets[:,:,i]=position_of_planets[:,:,i-1]+v*dt
			acceleration=self.compute_force_one_body(position_of_planets[:,:,i], self.mass_star)
			v+=acceleration*dt
			if i % int(N/float(time)) == 0 and i != 0:
				no_years+=1
				print "Done with year:", no_years
		times=np.linspace(0, time, N)
		self.mysolarsystem.orbitXml(position_of_planets, times)

		#Implement Euler-Cromer PLOT THIS TOO, COMPARE
		"""
		i=0
		no_years=0
		for i in xrange(1,N):
			acceleration=self.compute_force_one_body(position_of_planets[:,i-1,:], self.mass_star)
			v+=acceleration*dt
			position_of_planets[:,i,:]=position_of_planets[:,i-1,:]+v*dt
			if i % int(N/float(time)) == 0 and i != 0:
				no_years+=1
				print "Done with year:", no_years
		"""

		for k in xrange(self.no_planets):
			plt.plot(position_of_planets[0,k,:], position_of_planets[1,k,:])
			plt.hold('on')
		plt.plot(0,0, marker='*')
		plt.legend(xrange(self.no_planets))
		plt.xlabel('Distance from star (AU)')
		plt.ylabel('Distance from star (AU)')
		plt.grid('on')
		plt.axis('equal')
		plt.show()

		steps_per_year=int(N/float(time))

		self.mysolarsystem.checkPlanetPositions(position_of_planets, time, steps_per_year)


	def compute_orbit_most_massive_planets(self,time, N,no_massive_planets):
		if N%10 != 0:
			print "Need a timestep which is divisible by zero"
			sys.exit(1)

		masses_sorted=sorted(self.mass_of_planets, reverse=True)
		heaviest_planets=masses_sorted[:no_massive_planets]
		x0_heaviest=np.zeros(no_massive_planets+1)
		y0_heaviest=np.zeros(no_massive_planets+1)
		vx0_heaviest=np.zeros(no_massive_planets+1)
		vy0_heaviest=np.zeros(no_massive_planets+1)
		mass_heaviest=np.zeros(no_massive_planets+1)
		radius_heaviest=np.zeros(no_massive_planets+1)

		for m in xrange(len(heaviest_planets)):
			correct_planets_indices=((np.where(self.mass_of_planets==heaviest_planets[m])[0])[0])
			x0_heaviest[m]=self.x0[correct_planets_indices]
			y0_heaviest[m]=self.y0[correct_planets_indices]
			vx0_heaviest[m]=self.vx0[correct_planets_indices]
			vy0_heaviest[m]=self.vy0[correct_planets_indices]
			mass_heaviest[m]=self.mass_of_planets[correct_planets_indices]
			radius_heaviest[m]=self.radius_of_planets_in_AU[m]

		mass_heaviest[len(heaviest_planets)]=self.mass_star
		radius_heaviest[len(heaviest_planets)]=self.radius_of_star_in_AU
		no_planets=len(mass_heaviest)
		total_momentum_x=np.sum(vx0_heaviest*mass_heaviest)
		total_momentum_y=np.sum(vy0_heaviest*mass_heaviest)
		v0_sun_x=total_momentum_x/float(self.mass_star)
		v0_sun_y=total_momentum_y/float(self.mass_star)
		vx0_heaviest[-1]=-v0_sun_x
		vy0_heaviest[-1]=-v0_sun_y

		dt=time/float(N)
		position_of_planets=np.zeros(shape=(2, no_planets, int(N/10)))
		position_of_planets[0,:, 0]=x0_heaviest
		position_of_planets[1,:,0]=y0_heaviest

		v=np.zeros(shape=(2, no_planets))
		v[0,:]=vx0_heaviest
		v[1,:]=vy0_heaviest
		time_taken=np.linspace(0, time, N)



		no_years=0
		velocity_of_sun=np.zeros(shape=(int(N/10),2))
		velocity_of_sun[0,0]=vx0_heaviest[-1]
		velocity_of_sun[0,1]=vy0_heaviest[-1]
		temp_pos=np.zeros(shape=(2, no_planets))
		temp_pos[0,:]=x0_heaviest
		temp_pos[1,:]=y0_heaviest
		acceleration=self.compute_force_two_bodies(temp_pos, mass_heaviest)
		v+=0.5*acceleration*dt

		for i in xrange(1,N):
			temp_pos+=v*dt
			acceleration=self.compute_force_two_bodies(temp_pos, mass_heaviest)
			v+=acceleration*dt
			if i % 10 == 0:
				position_of_planets[:, :, int(i/10)]=temp_pos
				velocity_of_sun[int(i/10),:]=v[:, -1]
			if i % int(N/float(time)) == 0 and i != 0:
				no_years+=1
				print "Done with year:", no_years

		luminosity=np.ones(int(N/10))
		t=np.linspace(0, time, int(N/10))


		for l in range(no_planets-1):
			planet_in_front=np.logical_and(position_of_planets[0, l, :]> -self.radius_of_star_in_AU, position_of_planets[0, l, :] < self.radius_of_star_in_AU)
			luminosity[planet_in_front] = luminosity[planet_in_front]*(1-((radius_heaviest[l])**2/(float(radius_heaviest[-1])**2)))
		plt.plot(t, luminosity)
		plt.show()




		for k in xrange(no_planets):
			plt.plot(position_of_planets[0,k,:], position_of_planets[1,k,:])
			plt.hold('on')
		plt.grid('on')
		plt.axis('equal')
		plt.xlabel('Distance from star (AU)')
		plt.ylabel('Distance from star (AU)')
		plt.show()

		v_sun=np.linalg.norm(velocity_of_sun, axis=1)
		v_r_max=velocity_of_sun[:,0].max()
		random_noise=np.random.normal(0, v_r_max/5.0, int(N/10))
		v_r_with_random=velocity_of_sun[:,0]+random_noise
		#plt.plot(t, v_sun)
		#plt.hold("on")
		plt.plot(t, v_r_with_random)
		plt.hold("on")
		#plt.plot(t, velocity_of_sun[:,1])
		#plt.legend(["v", "v_x", "v_y"])
		plt.show()




	def compute_orbit_all_planets(self, time, N):
		dt=time/float(N)
		position_of_planets=np.zeros(shape=(self.no_planets+1, N, 2))
		position_of_planets[:-1,0, 0]=self.x0
		position_of_planets[:-1,0,1]=self.y0

		v=np.zeros(shape=(self.no_planets+1, 2))
		v[:-1,0]=self.vx0
		v[:-1,1]=self.vy0
		masses=np.zeros(len(self.mass_of_planets)+1)
		masses[:-1]=self.mass_of_planets
		masses[-1]=self.mass_star

		i=0
		no_years=0
		acceleration=self.compute_force_many_bodies(position_of_planets[:,0,:], masses)
		v+=0.5*acceleration*dt
		for i in xrange(1,N):
			position_of_planets[:,i,:]=position_of_planets[:,i-1,:]+v*dt
			acceleration=self.compute_force_two_bodies(position_of_planets[:,i,:], masses)
			v+=acceleration*dt
			if i % int(N/float(time)) == 0 and i != 0:
				no_years+=1
				print "Done with year:", no_years

		for k in xrange(self.no_planets+1):
			plt.plot(position_of_planets[k,:,0], position_of_planets[k,:,1])
			plt.hold('on')
		plt.grid('on')
		plt.xlabel('Distance from star (AU)')
		plt.ylabel('Distance from star (AU)')
		plt.show()
		position_of_planets=np.reshape(position_of_planets, (2, self.no_planets+1, N))















n_years=20
steps_per_year=20000
seed=6392
n_steps=steps_per_year*n_years
simulate_planets=Calculate_planet_orbits(seed)
simulate_planets.compute_orbits_one_body(n_years, n_steps, debug=False)
#simulate_planets.compute_orbit_most_massive_planets(n_years, n_steps, 3)
#simulate_planets.compute_orbit_all_planets(n_years, n_steps)
