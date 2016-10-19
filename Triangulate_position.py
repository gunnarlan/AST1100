import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.interpolate as inter
from AST1100SolarSystem import AST1100SolarSystem

class Launch_simulation:

    def __init__(self, no_years, steps_per_year, seed, target_planet=6):
		self.no_years=no_years
		self.steps_per_year=steps_per_year
		self.seed=seed
		self.target_planet=target_planet
		self.mysolarsystem=AST1100SolarSystem(self.seed,hasMoons=False)
		self.major_axis=self.mysolarsystem.a
		self.mass_star=self.mysolarsystem.starMass
		self.no_planets=self.mysolarsystem.numberOfPlanets
		self.mass_of_planets=self.mysolarsystem.mass
		self.radius_planets=self.mysolarsystem.radius
		self.period_planets=self.mysolarsystem.period
		self.vx0=self.mysolarsystem.vx0
		self.vy0=self.mysolarsystem.vy0
		[self.planetPos, t]=np.load('planetPositions.npy')
		self.positionFunction=inter.interp1d(t, self.planetPos)

    def get_position(self, t0):
        my_pos_x=np.random.uniform(-20, 20)
        my_pos_y=np.random.uniform(-20, 20)
        r=np.array([my_pos_x, my_pos_y])
        x_y=np.zeros(shape=(self.no_planets-1, 2))
        tol=1e-5
        diff=np.zeros(self.no_planets-1)
        for k in range(self.no_planets-1):
            r1=np.linalg.norm(r)
            r2=np.linalg.norm(r-self.positionFunction(t0)[:, k])
            r3=np.linalg.norm(r-self.positionFunction(t0)[:, k+1])
            x1=0
            y1=0
            x2=self.positionFunction(t0)[0,k]
            y2=self.positionFunction(t0)[1,k]
            x3=self.positionFunction(t0)[0,k+1]
            y3=self.positionFunction(t0)[1, k+1]
            x,y,difference=self.triangulate_analytic(x1,y1,r1,x2,y2,r2,x3,y3,r3)
            x_y[k, 0]=x
            x_y[k, 1]=y
            diff[k]=difference
        if (diff > tol).any():
            print diff.max()
            print "Oh no, one failed :("
            sys.exit(1)
        print "My pos x:", my_pos_x
        print "My pos y:", my_pos_y
        #return x1, y1, r1, x2, y2, r2, x3, y3, r3


    def check_position_from_array(self, t0):
	position_array=np.load('pos.npy')
	print position_array
	x_y=np.zeros(shape=(self.no_planets-1, 2))
        tol=1e-5
        diff=np.zeros(self.no_planets-1)
	for k in range(self.no_planets-1):
            r1=position_array[-1]
	    r2=position_array[k]
            r3=position_array[k+1]
            x1=0
            y1=0
            x2=self.positionFunction(t0)[0,k]
            y2=self.positionFunction(t0)[1,k]
            x3=self.positionFunction(t0)[0,k+1]
            y3=self.positionFunction(t0)[1, k+1]
            x,y,difference=self.triangulate_analytic(x1,y1,r1,x2,y2,r2,x3,y3,r3)
            x_y[k, 0]=x
            x_y[k, 1]=y
            diff[k]=difference
        if (diff > tol).any():
            print diff.max()
            print "Oh no, one failed :("
       	min_index=np.argmin(np.abs(diff))
	print "Lowest x:", x_y[min_index, 0]
	print "Lowest y:", x_y[min_index, 1]


    
    def triangulate_analytic_sun_at_center(self,r1,x2,y2,r2,x3,y3,r3):
        gamma=(r1**2+x2**2+y2**2-r2**2)/(2.0*x2)
	a=(y2**2)/(float(x2**2))
	b=-2.0*gamma*y2/x2
	c=gamma**2-r1**2
	y_plus=(-b+np.sqrt((b**2)-4*a*c))/(2.0*a)
	y_minus=(-b-np.sqrt((b**2)-4*a*c))/(2.0*a)
        x_plus=gamma-y_plus*y2/float(x2)
        x_minus=gamma-y_minus*y2/float(x2)
        difference_plus=(x_plus-x3)**2+(y_plus-y3)**2-r3**2
        difference_minus=(x_minus-x3)**2+(y_minus-y3)**2-r3**2
        if abs(difference_minus) < abs(difference_plus):
            print "Difference minus", difference_minus
            print x_minus, y_minus
            return x_minus, x_plus, difference_minus
        else:
            print "Difference plus", difference_plus
            print x_plus, y_plus
            return x_plus, y_plus, difference_plus



    def triangulate_analytic(self, x1,y1,r1,x2,y2,r2,x3,y3,r3):
        gamma=(r1**2+x2**2+y2**2-x1**2-y1**2-r2**2)/(float(2.0*(x2-x1)))-x1
        Lambda=(y1-y2)/float(x2-x1)
        y_plus=(-2.0*(gamma*Lambda-y1)+np.sqrt(4*(gamma*Lambda-y1)**2-4*(Lambda**2+1)*(gamma**2+y1**2-r1**2)))/float(2*(Lambda**2+1))
        y_minus=(-2.0*(gamma*Lambda-y1)-np.sqrt(4*(gamma*Lambda-y1)**2-4*(Lambda**2+1)*(gamma**2+y1**2-r1**2)))/float(2*(Lambda**2+1))
        x_plus=gamma+x1+(Lambda*y_plus)
        x_minus=gamma+x1+(Lambda*y_minus)
        difference_plus=(x_plus-x3)**2+(y_plus-y3)**2-r3**2
        difference_minus=(x_minus-x3)**2+(y_minus-y3)**2-r3**2
        if abs(difference_minus) < abs(difference_plus):
            print "Difference minus", difference_minus
            print x_minus, y_minus
            return x_minus, x_plus, difference_minus
        else:
            print "Difference plus", difference_plus
            print x_plus, y_plus
            return x_plus, y_plus, difference_plus


launch=Launch_simulation(100, 20000, 6392)
#launch.get_position(3.5909)

position_array=np.load('pos.npy')
#launch.triangulate_analytic_sun_at_center(position_array[-1],)
launch.check_position_from_array(3.601)
"""
for i in range(1000):
    launch.get_position(i/1000.0)
"""
