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
        self.mysystem.landOnPlanet(int(self.target_planet), 'Part7_Satellite_Instructions_landing.txt')








if __name__ == "__main__":
	launch_satellite=Real_launch_satellite(seed=6392)
	launch_satellite.launch()
