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
		self.AU=1.495979*1e11
		self.star_temperature=self.mysolarsystem.temperature
		self.radius_of_star_in_AU=1e3*self.mysolarsystem.starRadius*1.0/(self.AU)
		self.radius_of_planets_in_AU=1e3*self.mysolarsystem.radius*1.0/(self.AU)
		self.x0=self.mysolarsystem.x0
		self.y0=self.mysolarsystem.y0
		self.vx0=self.mysolarsystem.vx0
		self.vy0=self.mysolarsystem.vy0
		self.G=4*(np.pi)**2
		self.mu=float(self.mass_star*self.G)
		self.legend=['Byappo (0)', 'Domin (1)', 'Munnmon (2)', 'Pjeng (3)', 'Plaging (4)', 'Psiwe (5)', 'Hiffre (6)']

	def compute_force_one_body(self,position, m_object_producing_force):
		acceleration=np.zeros_like(position)
		for i in xrange(len(position[0])):
			distance=np.linalg.norm(position[:,i])
			a=-self.mu/(distance**3)
			a_vec=a*position[:,i]
			acceleration[:,i]=a_vec
		return acceleration

	def compute_habitable_r(self, T_min=260, T_max=390):
		r_min=(self.radius_of_star_in_AU*self.star_temperature**2)/float(2.0*T_max**2)
		r_max=(self.radius_of_star_in_AU*self.star_temperature**2)/float(2.0*T_min**2)
		return r_min, r_max

	def compute_force_one_body_with_energy(self,position, m_object_producing_force):
		acceleration=np.zeros_like(position)
		for i in xrange(len(position[0])):
			distance=np.linalg.norm(position[:,i])
			a=-self.mu/(distance**3)
			a_vec=a*position[:,i]
			acceleration[:,i]=a_vec
		return acceleration

	def compute_force_two_bodies(self, position, m_planets):
		acceleration=np.zeros_like(position)
		for i in xrange(len(m_planets)):
			a=np.zeros(2)
			if i != len(m_planets)-1:
				distance=position[:,i]-position[:,-1]
				a=-(self.mu*distance/(np.linalg.norm(distance)**3))
			else:
				for j in xrange(len(m_planets)-1):
					distance=position[:,j]-position[:,-1]
					a=a+(self.G*m_planets[j]/(float(np.linalg.norm(distance)**3)))*distance
			acceleration[:, i]=a
		return acceleration


	def compute_force_many_bodies(self, position, m_planets):
		acceleration=np.zeros_like(position)
		for i in xrange(len(m_planets)):
			a=np.zeros(2)
			for k in xrange(len(m_planets)):
				if k != i:
					distance=position[k]-position[i]
					a+=(self.G*m_planets[k]/(float(np.linalg.norm(distance))**3))*distance
			acceleration[i]=a
		return acceleration



	def  compute_orbits_one_body(self, time, N, debug=False, tol=0.01, upper=1, lower=1e-15, check_energy=False, plot_habitable=False):
		"""
		Ignore interplanetary forces and all forces on the star
		"""

		dt=time/float(N)
		position_of_planets=np.zeros(shape=(2, self.no_planets, N))
		position_of_planets[0,:,0]=self.x0
		position_of_planets[1,:,0]=self.y0

		if check_energy==True:
			energy=np.zeros(shape=(3, self.no_planets, N))
			x_vec=np.array([self.x0, self.y0])
			v_vec=np.array([self.vx0, self.vy0])
			energy[0,:,0]=-self.mass_of_planets*self.mu/(np.linalg.norm(x_vec, axis=0)) #Potential energy
			energy[1,:,0]=0.5*self.mass_of_planets*(np.linalg.norm(v_vec, axis=0))**2 #Kinetic energy
			energy[2,:,0]=np.cross(self.mass_of_planets*v_vec, x_vec, axis=0)

		v=np.zeros(shape=(2,self.no_planets))
		v[0,:]=self.vx0
		v[1,:]=self.vy0

		if debug:  #Check to find required h
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

		#Euler - to check leapfrog
		"""
		i=0
		no_years=0
		for i in xrange(1,N):
			acceleration=self.compute_force_one_body(position_of_planets[:,:,i-1], self.mass_star)
			v+=acceleration*dt
			position_of_planets[:,:,i]=position_of_planets[:,:,i-1]+v*dt
			energy[0,:,i]=-self.mass_of_planets*self.mu/(np.linalg.norm(position_of_planets[:,:,i], axis=0)) #Potential energy
			energy[1,:,i]=0.5*self.mass_of_planets*(np.linalg.norm(v, axis=0))**2 #Kinetic energy
			momentum=self.mass_of_planets*v
			energy[2, :, i]=np.cross(momentum, position_of_planets[:,:,i], axis=0)
			if i % int(N/float(time)) == 0 and i != 0:
				no_years+=1
				print "Done with year:", no_years

		"""
		#Leapfrog
		i=0
		no_years=0
		acceleration=self.compute_force_one_body(position_of_planets[:,:,0], self.mass_star)
		
		v+=0.5*acceleration*dt
		if check_energy==True:
			for i in xrange(1,N):
				position_of_planets[:,:,i]=position_of_planets[:,:,i-1]+v*dt
				acceleration=self.compute_force_one_body(position_of_planets[:,:,i], self.mass_star)
				v+=acceleration*dt
				energy[0,:,i]=-self.mass_of_planets*self.mu/(np.linalg.norm(position_of_planets[:,:,i], axis=0)) #Potential energy
				energy[1,:,i]=0.5*self.mass_of_planets*(np.linalg.norm(v, axis=0))**2 #Kinetic energy
				momentum=self.mass_of_planets*v
				energy[2, :, i]=np.cross(momentum, position_of_planets[:,:,i], axis=0)
				
				if i % int(N/float(time)) == 0 and i != 0:
					no_years+=1
					print "Done with year:", no_years

		else:
			for i in xrange(1,N):
				position_of_planets[:,:,i]=position_of_planets[:,:,i-1]+v*dt
				acceleration=self.compute_force_one_body(position_of_planets[:,:,i], self.mass_star)
				v+=acceleration*dt
				if i % int(N/float(time)) == 0 and i != 0:
					no_years+=1
					print "Done with year:", no_years

		times=np.linspace(0, time, N)
		#self.mysolarsystem.orbitXml(position_of_planets, times)

		if check_energy==True:
			total_energy=energy[0,:,:]+energy[1,:,:]
			for l in xrange(self.no_planets):
				relative_error_energy=100*(total_energy[l,:]-total_energy[l,0])/(float(total_energy[l,0]))
				print "Energy: maximum error of planet %d in percent is %.6e" %(l, max(abs(relative_error_energy)))
				plt.plot(times,relative_error_energy)
				plt.hold('on')
			plt.xlabel('Time [years]')
			plt.ylabel('Relative error of total energy [%]')
			plt.grid('on')
			plt.title('Relative error in total energy as a function of time, with dt='+str(dt))
			plt.legend(xrange(self.no_planets))
			plt.show()
			plt.clf()

			for m in xrange(self.no_planets):
				relative_error_angular_momentum=100*(energy[2,m,:]-energy[2,m,0])/(float(energy[2,m,0]))
				print "Angular momentum: maximum error of energy planet %d in percent is %.4e" %(m, max(abs(relative_error_angular_momentum)))
				plt.plot(times,relative_error_angular_momentum)
				plt.hold('on')
			plt.xlabel('Time [years]')
			plt.ylabel('Relative error of angular momentum [%]')
			plt.grid('on')
			plt.title('Relative error of angular momentum as a function of time, with dt='+str(dt))
			plt.legend(xrange(self.no_planets))
			plt.show()

		if plot_habitable==True:
			r_min, r_max=self.compute_habitable_r()
			no_of_r=500
			theta=np.linspace(0, 2.0*np.pi, 100000)
			r=np.linspace(r_min, r_max, no_of_r)
			for j in range(no_of_r):
				plt.plot(r[j]*np.sin(theta), r[j]*np.cos(theta), 'r', label='_nolegend_')
				plt.hold('on')
			print "Maximum habital r [AU]:", r_max
			print "Minimum habitable r [AU]:", r_min
			plt.title('Habitable zone')
		else:
			plt.title('Orbit of the planets in our solar system')

		for k in xrange(self.no_planets):
			plt.plot(position_of_planets[0,k,:], position_of_planets[1,k,:])
			plt.hold('on')
		plt.plot(0,0, marker='*')
		plt.legend(self.legend, prop={'size':10})
		plt.xlabel('Distance from star [AU]')
		plt.ylabel('Distance from star [AU]')
		plt.grid('on')
		plt.axis('equal')
		plt.show()

		steps_per_year=int(N/float(time))

		#self.mysolarsystem.checkPlanetPositions(position_of_planets, time, steps_per_year)


	def compute_orbit_most_massive_planets(self,time, N,no_massive_planets):
		"""
		Computes the radial velocity by taking into account the most massive planets
		"""

		if N%10 != 0:
			print "Need a timestep which is divisible by ten"
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

		radius_heaviest[len(heaviest_planets)]=self.radius_of_star_in_AU
		no_planets=len(mass_heaviest)
		total_momentum_x=np.sum(vx0_heaviest*mass_heaviest)
		total_momentum_y=np.sum(vy0_heaviest*mass_heaviest)
		v0_sun_x=total_momentum_x/float(self.mass_star)
		v0_sun_y=total_momentum_y/float(self.mass_star)
		vx0_heaviest[-1]=-v0_sun_x
		vy0_heaviest[-1]=-v0_sun_y

		mass_heaviest[len(heaviest_planets)]=self.mass_star
		

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


		for k in xrange(no_planets):
			plt.plot(position_of_planets[0,k,:], position_of_planets[1,k,:])
			plt.hold('on')
		plt.grid('on')
		plt.title('Motion around the center of mass in our solar system')
		plt.axis('equal')
		plt.grid('on')
		plt.xlabel('Distance from initial position of the star [AU]')
		plt.ylabel('Distance from initial position of the star [AU]')
		plt.show()

		t=np.linspace(0, time, len(v_sun))
		v_r_max=velocity_of_sun[:,0].max()
		random_noise=np.random.normal(0, v_r_max/5.0, int(N/10))
		v_r_with_random=velocity_of_sun[:,0]+random_noise
		#plt.plot(t, v_sun)
		#plt.hold("on")
		plt.plot(t, v_r_with_random)
		plt.title('Radial speed of our star around the center of mass, with noise')
		plt.xlabel('Time [year]')
		plt.ylabel('Speed of star [AU/year]')
		plt.grid('on')
		plt.show()




	def compute_orbit_all_planets(self, time, N):
		"""
		Unused function implemented as reference
		"""

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
		plt.xlabel('Distance from star [AU]')
		plt.ylabel('Distance from star [AU]')
		plt.show()
		position_of_planets=np.reshape(position_of_planets, (2, self.no_planets+1, N))














if __name__ == "__main__":
	n_years=100
	steps_per_year=20000
	seed=6392
	n_steps=steps_per_year*n_years
	simulate_planets=Calculate_planet_orbits(seed)
	simulate_planets.compute_orbits_one_body(n_years, n_steps, debug=False)
	#simulate_planets.compute_orbit_most_massive_planets(n_years, n_steps, 3)

