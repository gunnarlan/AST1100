import numpy as np
import matplotlib.pyplot as plt
import sys
from PIL import Image
from AST1100SolarSystem import AST1100SolarSystem


class Orientate_satellite:
    def __init__(self):
        self.AU=1.496e11
        self.theta_0=np.pi/2.0
        self.infile=open('himmelkule.npy', 'rb')
        self.celestial_sphere=np.load(self.infile)
        self.infile.close()

    def find_pixels(self):
	"""
	Find size of image
	"""
        ref_image=Image.open('sample0000.png')
        imarray=np.array(ref_image)
        ref_image.close()
        self.number_of_pix=imarray.shape
        print self.number_of_pix
        ref_image=None
        imarray=None

    def compute_extreme_values(self, alpha_phi, alpha_theta):
	"""
	Finds the extrme values of the tangent coordinates
	"""
        extreme_values=np.zeros(4)
        numerator_x=2*np.sin(alpha_phi/2.0)
        denominator_x=1+np.cos(alpha_phi/2.0)
        numerator_y=2*np.sin(alpha_theta/2.0)
        denominator_y=1+np.cos(alpha_theta/2.0)
        extreme_values[0]=-numerator_x/float(denominator_x)
        extreme_values[1]=numerator_x/float(denominator_x)
        extreme_values[2]=numerator_y/float(denominator_y)
        extreme_values[3]=-numerator_y/float(denominator_y) #Flipping the y-axis
        return extreme_values

    def compute_theta_phi(self, x, y, phi_0, theta_0=np.pi/2):
	"""
	Computes the angles associated with the reference sphere
	"""
        angles=np.zeros(shape=[x.shape[0], y.shape[0],2])
        xx, yy=np.meshgrid(x,y)
        rho=np.sqrt(xx**2+yy**2)
        c=2.0*np.arctan(rho/2.0)
        theta=theta_0-np.arcsin(np.cos(c)*np.cos(theta_0)+yy*np.sin(c)*np.sin(theta_0)/rho)
        phi=phi_0+np.arctan(xx*np.sin(c)/(rho*np.sin(theta_0)*np.cos(c)-yy*np.cos(theta_0)*np.sin(c)))
        angles[:,:,0]=theta
        angles[:,:,1]=phi
        return angles


    def compute_reference_sphere(self, x, y):
	"""
	Computes and saves the reference sphere
	"""
        theta_0=np.pi/2.0
        pictures=np.zeros(shape=(360, int(y.shape[0]), int(x.shape[0]), 3), dtype=np.uint8)
        xx, yy=np.meshgrid(x,y)
        rho=np.sqrt(xx**2+yy**2)
        c=2.0*np.arctan(rho/2.0)
        theta=theta_0-np.arcsin(yy*np.sin(c)/rho)
        for phi_0 in range(0, 360):
            phi_in_rad=np.deg2rad(phi_0)
            phi=phi_in_rad+np.arctan(xx*np.sin(c)/(rho*np.cos(c)))
            for k in range(len(x)):
                for j in range(len(y)):
                    pixnum=AST1100SolarSystem.ang2pix(theta[j,k], phi[j,k])
                    temp=self.celestial_sphere[pixnum]
                    pictures[phi_0, j, k, :]=[temp[2], temp[3], temp[4]]
            print "Done with phi: ", phi_0
        np.save("Reference_sphere.npy", pictures)

    def compute_projections(self, alpha_phi, alpha_theta):
	"""
	Sets up the reference sphere computations
	"""
        self.find_pixels()
        extreme_values=self.compute_extreme_values(alpha_phi, alpha_theta)
        x=np.linspace(extreme_values[0], extreme_values[1], self.number_of_pix[1])
        y=np.linspace(extreme_values[2], extreme_values[3], self.number_of_pix[0])
        self.compute_reference_sphere(x,y)

    def compute_phi(self, input_image):
	"""
	Uses least square to find the best phi
	"""
        ref_image=Image.open(input_image)
        imarray=np.array(ref_image)
        ypix=imarray[0,:,0].shape
        ref_image.close()

        infile=open('Reference_sphere.npy', 'rb')
        reference_sphere=np.load(infile)
        infile.close()
        print reference_sphere.shape

        diff = 1000000000
        j=0
        for k in range(reference_sphere.shape[0]):
            diff_img=(imarray-reference_sphere[k, :, :, :])**2
            least_square=np.sum(diff_img)
            print least_square
            if least_square < diff:
                diff = least_square
                j=k
        img3=Image.fromarray(reference_sphere[j, :, :, :])
        img3.save('Compute_image.png')
        return j, diff

    def check_angle(self):
	"""
	Warpper for compute_phi function
	"""
        self.find_pixels()
        alpha_theta=np.deg2rad(70)
        alpha_phi=np.deg2rad(70)
        extreme_values=self.compute_extreme_values(alpha_phi, alpha_theta)
        x=np.linspace(extreme_values[0], extreme_values[1], self.number_of_pix[1])
        y=np.linspace(extreme_values[2], extreme_values[3], self.number_of_pix[0])
        phi_0=20
        phi_0=np.deg2rad(phi_0)
        j, diff=self.compute_phi("find_orient.png")
        print "j=", j
        print "diff=", diff




if __name__ == "__main__":
	orient=Orientate_satellite()
	alpha_phi_radians=np.deg2rad(70)
	#orient.compute_projections(alpha_phi_radians, alpha_phi_radians)
	orient.check_angle()
	#j, diff=orient.compute_phi('sample0000.png')
	#print "j=", j
	#print "diff=", diff
