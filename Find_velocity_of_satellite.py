import numpy as np
import matplotlib.pyplot as plt
import sys
from AST1100SolarSystem import AST1100SolarSystem

class Find_Velocity():

    def __init__(self):
        self.solar_system=AST1100SolarSystem(6392)
        self.c=299792458
        self.solar_system.getRefStars()
        self.m_s_to_AU_year=(365*24*60*60)/(1.49597871*1e11)

    def compute_radial_velocity(self,delta_lambda_1, delta_lambda_2, H_alpha=656.3*1e-9):
        radial_velocity_1=self.c*(delta_lambda_1)/(float(H_alpha))
        radial_velocity_2=self.c*(delta_lambda_2)/(float(H_alpha))
        return radial_velocity_1, radial_velocity_2

    def compute_velocity_of_sat(self, d_lamb_wrt_star_1, d_lamb_wrt_star_2, d_lamb_wrt_sat_1, d_lamb_wrt_sat_2, phi_1, phi_2):
        if abs(phi_1-phi_2)<1e-18:
            print "Warning, your stars are collinear. There is not enough information to determine the velocity"
            sys.exit(1)
        v_refstar_1, v_refstar_2=self.compute_radial_velocity(d_lamb_wrt_star_1, d_lamb_wrt_star_2)
        print "Velocity of star 1:",  v_refstar_1
        print "Velocity of star 2:",  v_refstar_2
        v_rel_1, v_rel_2=self.compute_radial_velocity(d_lamb_wrt_sat_1, d_lamb_wrt_sat_2)
        v_sat=np.array(([v_refstar_1-v_rel_1, v_refstar_2-v_rel_2]))
        prefac=1/float(np.sin(phi_2-phi_1))
        A_inv=np.array(([np.sin(phi_2), -np.sin(phi_1)], [-np.cos(phi_2), np.cos(phi_1)]))
        v=self.m_s_to_AU_year*prefac*(A_inv.dot(v_sat))
        print v



lambda_1_wrt_star=-0.020159867169*1e-9
lambda_2_wrt_star=-0.014180672906*1e-9
lambda_1_wrt_sat=-0.041737905225*1e-9
lambda_2_wrt_sat=-0.082521637840*1e-9
phi_1=np.deg2rad(171.006135)
phi_2=np.deg2rad(87.754469)
velocity=Find_Velocity()
velocity.compute_velocity_of_sat(lambda_1_wrt_star, lambda_2_wrt_star, lambda_1_wrt_sat, lambda_2_wrt_sat, phi_1, phi_2)
