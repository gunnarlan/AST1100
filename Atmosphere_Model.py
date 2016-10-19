import numpy as np
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem
from scipy.optimize import brentq

class Planet_Atmosphere():
	def __init__(self, seed, mu, gamma):
		self.gamma=gamma
		self.mu=mu
		self.mysolarsystem=AST1100SolarSystem(seed)
		self.AU=1.495978707*1e11
		self.G=6.674*1e-11
		self.target_planet=6
		self.M_sun_in_kg=1.989*1e30
		self.k= 1.380648*1e-23
		self.radius_planet_in_m=1e3*self.mysolarsystem.radius[self.target_planet]
		print self.radius_planet_in_m/float(1e3)
		self.M_in_kg=self.M_sun_in_kg*(self.mysolarsystem.mass[self.target_planet])
		self.rho_0=self.mysolarsystem.rho0[self.target_planet]
		self.atomic_mass_unit=1.661*1e-27
		self.T_0=(self.mysolarsystem.temperature)*np.sqrt(self.mysolarsystem.starRadius/(2.0*self.mysolarsystem.a[self.target_planet]*self.AU/(1000.0)))
		self.mass_satellite=1100
		self.omega=(2.0*np.pi)/(24*60*60*self.mysolarsystem.period[self.target_planet])
		self.A=6

	def compute_r_T_2(self):
		r_inv=(1.0/self.radius_planet_in_m)-(self.T_0/2.0)*(self.gamma/(self.gamma-1))*(self.k/(self.G*self.M_in_kg*self.mu*self.atomic_mass_unit))
		return (r_inv)**(-1)

	def compute_atmosphere_adiabatic_analytic(self, r):
		eps=1e-8
		if isinstance(r, (list, tuple, np.ndarray)):
			if r[0] < self.radius_planet_in_m-eps:
				print "Uh-oh, you choose an r inside the planet"
				sys.exit(1)
		else:
			if r < self.radius_planet_in_m-eps:
				print "Uh-oh, you choose an r inside the planet"
				sys.exit(1)

		prefac=((self.gamma-1.0)/(self.gamma))*((self.G*self.M_in_kg*self.mu*self.atomic_mass_unit)/(self.k))
		T=prefac*((1.0/r)-(1.0/self.radius_planet_in_m))+self.T_0
		rho=self.rho_0*((T/float(self.T_0))**(1.0/(self.gamma-1)))
		return rho

	def compute_atmosphere_istohermal_analytic(self,r_0, r):
		eps=1e-5
		if isinstance(r, (list, tuple, np.ndarray)):
			if r[0] < r_0-eps:
				print "Uh-oh, tried to calculate isothermal atmosphere in adiabatic region"
				sys.exit(1)
		else:
			if r < r_0-eps:
				print "Uh-oh, tried to calculate isothermal atmosphere in adiabatic region"
				sys.exit(1)

		rho_0=self.compute_atmosphere_adiabatic_analytic(r_0)
		exponent=((2.0*self.G*self.M_in_kg*self.mu*self.atomic_mass_unit)/(self.k*self.T_0))
		rho=rho_0*np.exp(exponent*(1.0/r-1.0/r_0))
		return rho

	def compute_critical_r(self, eps):
		treshold=1e-15
		r_isothermal=self.compute_r_T_2()
		rho_0=self.compute_atmosphere_adiabatic_analytic(r_isothermal)
		rho_fac=np.log((eps/float(rho_0)))
		prefac=self.T_0*self.k/(2.0*self.G*self.M_in_kg*self.mu*self.atomic_mass_unit)
		r=((prefac*rho_fac)+(1.0/r_isothermal))**(-1)
		if self.compute_atmosphere_istohermal_analytic(r_isothermal, r) > eps+treshold:
			print "Uh-oh, critical radius wrongly impleneted, exiting"
			sys.exit(1)
		else:
			return r

	def compute_closest_r(self, r):
		r_isothermal=self.compute_r_T_2()
		rho_0=self.compute_atmosphere_adiabatic_analytic(r_isothermal)
		factor_1=self.G*self.M_in_kg/(500.0*self.A)
		exponent=2.0*self.G*self.M_in_kg*self.mass_satellite*self.mu*self.atomic_mass_unit/(self.k*self.T_0)
		exponent=exponent*(1.0/r-1.0/r_isothermal)
		factor_2=rho_0*(np.exp(exponent))*(r**2)*(np.sqrt(float(self.G)*self.M_in_kg/r)-self.omega*r)**2
		return factor_1-factor_2

	def find_root_closest_r(self):
		r_isothermal=self.compute_r_T_2()
		print r_isothermal/(float(1e3))
		x_0=brentq(self.compute_closest_r, a=r_isothermal, b=1e8)
		return x_0


	def compute_atmosphere_analytic(self, r_isothermal, r):
		if r < r_isothermal:
			return self.compute_atmosphere_adiabatic_analytic(r)
		else:
			return self.compute_atmosphere_istohermal_analytic(r_isothermal,r)


	def compute_atmosphere_numeric(self):
		P_0=(self.k*self.rho_0*self.T_0)/(self.mu*self.atomic_mass_unit)
		adiabtic_constant=(P_0**(1-self.gamma))*(self.T_0**(self.gamma))
		r=self.radius_planet_in_m

		T=self.T_0
		P=P_0
		rho=[self.rho_0]
		ideal_gas_constant=(self.mu*self.atomic_mass_unit)/float(self.k)

		exponent=(self.gamma-1)/(float(self.gamma))
		mu_gravity=-self.G*self.M_in_kg

		dr=1
		while T > self.T_0/2.0:
			P+=mu_gravity*rho[-1]/(float(r)**2)*dr
			T=(float(adiabtic_constant)/(P**(1-self.gamma)))**(1.0/self.gamma)
			rho.append(ideal_gas_constant*(P/float(T)))
			r+=dr

		no_points_adibatic=len(rho)
		max_r_adiabatic=r
		print max_r_adiabatic
		max_r_analytic=self.compute_r_T_2()
		print max_r_analytic


		while rho[-1]>1e-10:
			P+=mu_gravity*rho[-1]/(float(r)**2)*dr
			rho.append(ideal_gas_constant*(P/float(T)))
			r+=dr

		r_tot=np.linspace(self.radius_planet_in_m, r, len(rho))

		no_points_iso=len(rho)-no_points_adibatic
		r_iso=np.linspace(max_r_adiabatic, r,no_points_iso)
		r_adiabatic=np.linspace(self.radius_planet_in_m, max_r_adiabatic, no_points_adibatic)

		rho_adibatic_analytic=self.compute_atmosphere_adiabatic_analytic(r_adiabatic)
		rho_isothermal_analytic=self.compute_atmosphere_istohermal_analytic(max_r_adiabatic, r_iso)
		rho_tot_analytic=np.append(rho_adibatic_analytic, rho_isothermal_analytic)
		rho_diff=rho_tot_analytic-rho
		print max(rho_diff)

		plt.plot(r_tot/(float(self.radius_planet_in_m)), rho)
		plt.hold('on')
		plt.plot(r_adiabatic/(float(self.radius_planet_in_m)), rho_adibatic_analytic)
		plt.plot(r_iso/(float(self.radius_planet_in_m)), rho_isothermal_analytic)
		plt.grid('on')
		plt.show()






if __name__=='__main__':
	masses_present=[28.0101]
	mu=np.mean(masses_present)
	print mu
	atmosphere=Planet_Atmosphere(seed=6392, mu=mu, gamma=1.4)
	#rho=atmosphere.compute_atmosphere_adiabatic_analytic(2057857)
	#rho=atmosphere.compute_atmosphere_istohermal_analytic()
	#rho=atmosphere.compute_atmosphere_numeric()
	r=atmosphere.compute_critical_r(1e-10)
	atmosphere.find_root_closest_r()

