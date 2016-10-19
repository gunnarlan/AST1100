import numpy as np
import scipy.interpolate as inter
from AST1100SolarSystem import AST1100SolarSystem

class Real_launch_satellite():

    def __init__(self, seed):
        self.mysystem=AST1100SolarSystem(seed)
	self.mass_of_planets=self.mysystem.mass
	[self.planetPos, t]=np.load('planetPositions.npy')
	self.positionFunction=inter.interp1d(t, self.planetPos)
	self.target_planet=6
	self.G=4*(np.pi)**2


    def compute_velocity_planets(self, t, dt):
	v=(self.positionFunction(t+dt)-self.positionFunction(t-dt))/float(2*dt)
	return v



    def start_orbital_injection(self, t, position, velocity):
	planet_vel=self.compute_velocity_planets(t, 1e-8)
	relative_velocity=velocity-planet_vel[:, self.target_planet]
	position_vector=position-self.positionFunction(t)[:, self.target_planet]
	v_so=np.sqrt(self.mass_of_planets[self.target_planet]*self.G/float(np.linalg.norm(position_vector)))
	angle=np.arctan2(position_vector[1], position_vector[0])
	injection_velocity_x=v_so*np.sin(angle)-relative_velocity[0]
	injection_velocity_y=-v_so*np.cos(angle)-relative_velocity[1]
	injection_velocity=np.array([injection_velocity_x, injection_velocity_y])
	print "Injection velocity:", injection_velocity
	return injection_velocity
		

    def launch(self):
        self.mysystem.sendSatellite('Satellite_instructions.txt')
	inject_vel=self.start_orbital_injection(7.12474, [-7.41885809782,  -0.23167848846], [-0.241675952253 ,  -2.88286430256])
	t1=7.09589438124
	t2=7.5
	t3=8.5
	t4=9.0
	distance_1=np.linalg.norm([-7.4144299907 , -0.149121699967]-self.positionFunction(t1)[:, self.target_planet])  #0.0241044307971
	distance_2=np.linalg.norm([-7.27986774832,  -1.60803546091]-self.positionFunction(t2)[:, self.target_planet]) #7.12473 0.000425824411789
	distance_3=np.linalg.norm([-5.73753843321 , -4.88813682393]-self.positionFunction(t3)[:, self.target_planet])
	distance_4=np.linalg.norm([-4.42891695786 , -6.14234439191]-self.positionFunction(t4)[:, self.target_planet])
	print "Distance at time "+str(t1)+" equals: "+str(distance_1)
	print "Distance at time "+str(t2)+" equals: "+str(distance_2)
	print "Distance at time "+str(t3)+" equals: "+str(distance_3)
	print "Distance at time "+str(t4)+" equals: "+str(distance_4)












launch_satellite=Real_launch_satellite(seed=6392)
launch_satellite.launch()


#boost 7.09625056728 0.07187354 -0.87528548
