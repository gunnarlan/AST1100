import numpy as np
import matplotlib.pyplot as plt
import sys

class Find_Atmosphere():

	def __init__(self, seed, lambda_0, mass_of_molecule):
		self.seed=6392
		self.c=299792458
		self.lambda_0=lambda_0
		self.mass_of_molecule=mass_of_molecule
		boundary_lambda=(0.1/600.0)*self.lambda_0
		if self.lambda_0 < 600+boundary_lambda:
			print ("Error, need at least %.1f lambdas above lower boundary" %boundary_lambda)
			sys.exit(1)
		self.spectrum=np.load('spectrum_seed_92.npy')
		self.lower_index_of_interest=np.argmax(self.spectrum[:, 0]>lambda_0-boundary_lambda)
		self.upper_index_of_interest=np.argmax(self.spectrum[:, 0]>lambda_0+boundary_lambda)
		self.noise=np.load('sigma_noise.npy')


	def compute_error_in_sigma(self, T_upper=450, T_lower=150):
		avogadro=6.0221409*1e23
		k= 1.380648*1e-23
		def comp_FWHM(T):
			prefac=(2.0*self.lambda_0/self.c)
			root=(2.0*k*T*np.log(2))/(self.mass_of_molecule/(1000.0*avogadro))
			return prefac*np.sqrt(root)
		return comp_FWHM(T_lower), comp_FWHM(T_upper)


	def compute_error_in_lambda_center(self, max_velocity):
		return (max_velocity/float(self.c))*self.lambda_0

	def compute_error_in_F_min(self):
		return 0.3

	def Chi_square(self):
		F_min_min=0.7
		lower_sigma_error, upper_sigma_error=self.compute_error_in_sigma()
		lambda_error=self.compute_error_in_lambda_center(10000)
		F_error=self.compute_error_in_F_min()
		lambda_center=np.linspace(self.lambda_0-lambda_error, self.lambda_0+lambda_error, 300)
		F_min=np.linspace(F_min_min, F_min_min+F_error, 30)
		sigma=np.linspace(lower_sigma_error/20.0, 2.0*upper_sigma_error, 30)
		interesting_lambdas=self.spectrum[self.lower_index_of_interest:self.upper_index_of_interest,0]
		interesting_flux=self.spectrum[self.lower_index_of_interest:self.upper_index_of_interest,1]
		relevant_noise=self.noise[self.lower_index_of_interest:self.upper_index_of_interest, 1]

		chi_square=np.zeros(shape=(len(F_min), len(sigma),len(lambda_center)))
		divident_exp=2.0*(sigma)**2
		prefac=F_min-1
		for i in range(len(F_min)):
			for k in range(len(lambda_center)):
				temp=1+prefac[i]*np.exp(-(interesting_lambdas[None, :]-lambda_center[k])**2/(divident_exp[:, None]))
				chi_square[i,:,k]=np.sum(((interesting_flux[None, :]-temp[:,:])**2/(((relevant_noise[None, :])**2))), axis=1)
			print i


		smallest_chi=1e10
		index_F=1000
		index_sigma=1000
		index_lambda=1000
		for i in range(len(F_min)):
			for j in range(len(sigma)):
				for k in range(len(lambda_center)):
					if chi_square[i,j,k] < smallest_chi:
						smallest_chi=chi_square[i,j,k]
						index_F=i
						index_sigma=j
						index_lambda=k

		outfile=open('Chi_square_parameters.txt', 'a')
		outfile.write("%.5f  %.5f    %.5f   %.5f   %.5f   %.5f\n" %(self.lambda_0, chi_square[index_F, index_sigma, index_lambda], F_min[index_F], sigma[index_sigma], lambda_center[index_lambda], 				self.mass_of_molecule))
		outfile.close()

if __name__ == "__main__":
	mass_O2= 31.9988
	mass_H2O=18.0153
	mass_CO2=44.0095
	mass_CH4=16.0425
	mass_CO=28.0101
	mass_N2O=44.0128

	lambda_O2=[630, 690, 760]
	lambda_H2O=[720, 820, 940]
	lambda_CO2=[1400, 1600]
	lambda_CH4=[1660, 2200]
	lambda_CO=[2340]
	lambda_N2O=[2870]


	for k in range(len(lambda_O2)):
		compute_atomsphere=Find_Atmosphere(6392, lambda_O2[k], mass_O2)
		compute_atomsphere.Chi_square()

	for k in range(len(lambda_H2O)):
		compute_atomsphere=Find_Atmosphere(6392, lambda_H2O[k], mass_H2O)
		compute_atomsphere.Chi_square()

	for k in range(len(lambda_CO2)):
		compute_atomsphere=Find_Atmosphere(6392, lambda_CO2[k], mass_CO2)
		compute_atomsphere.Chi_square()

	for k in range(len(lambda_CH4)):
		compute_atomsphere=Find_Atmosphere(6392, lambda_CH4[k], mass_CH4)
		compute_atomsphere.Chi_square()

	for k in range(len(lambda_CO)):
		compute_atomsphere=Find_Atmosphere(6392, lambda_CO[k], mass_CO)
		compute_atomsphere.Chi_square()

	for k in range(len(lambda_N2O)):
		compute_atomsphere=Find_Atmosphere(6392, lambda_N2O[k], mass_N2O)
		compute_atomsphere.Chi_square()
