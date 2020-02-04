import numpy as np
import scipy.interpolate as inter
from AST1100SolarSystem import AST1100SolarSystem

class Real_launch_satellite():

    def __init__(self, seed):
        self.mysystem=AST1100SolarSystem(seed)
        self.mass_of_planets=self.mysystem.mass
        self.target_planet=6
	self.radius_planet_in_m=1e3*self.mysystem.radius[self.target_planet]
        self.G=6.674*1e-11




    def launch(self):
        self.mysystem.landOnPlanet(int(self.target_planet), 'Part6_Satellite_instructions_landing.txt')
	print (np.linalg.norm(np.array([ -1305086.79918075, -1971490.62974993,   729763.1460648 ]))-self.radius_planet_in_m)/1000.0
	print (np.linalg.norm(np.array([-1288663.74547901, -1946677.741875,    -819435.27243301]))-self.radius_planet_in_m)/1000.0
	print (np.linalg.norm(np.array([ 1082015.59827507,  1634509.40019873,  1510164.38669593 ]))-self.radius_planet_in_m)/1000.0
	print (np.linalg.norm(np.array([856437.88401898,  1293756.12314663, -1927094.20118337]))-self.radius_planet_in_m)/1000.0
	print (np.linalg.norm(np.array([ -900801.20857297, -1360772.0316381,   1860109.22459017 ]))-self.radius_planet_in_m)/1000.0
	




if __name__ == "__main__":
	launch_satellite=Real_launch_satellite(seed=6392)
	launch_satellite.launch()
