import numpy as np
import matplotlib.pyplot as plt

def compute_velocity_of_satellite(delta_lambda, lambda_0):
	prefac=delta_lambda/float(lambda_0)
	c=299792458
	return prefac*c

def compute_temperature_of_planet(sigma, lambda_0, m):
	c=299792458
	k=1.380648*1e-23
	atomic_mass_unit=1.661*1e-27
	mass=m*atomic_mass_unit
	FWHM=sigma*np.sqrt(8*np.log(2))
	T=(mass/(2.0*k*np.log(2)))*((FWHM*c/(2.0*lambda_0))**2)
	return T


infile=open("Chi_square_parameters.txt", "r")
infile_read=infile.readlines()
spectrum=np.load('spectrum_seed_92.npy')
no_wavelengths=np.zeros(len(infile_read)-1)
masses=np.zeros(len(infile_read)-1)
k=0
species=['O2','O2','O2','H2O','H2O','H2O','CO2','CO2','CH4', 'CH4', 'CO', 'N2O' ]
for line in infile_read[:-1]:
	lines=line.split()
	lambda_0=float(lines[0])
	lambda_estimate=lines[4]
	delta_lambda=float(lambda_estimate)-lambda_0
	velocity=compute_velocity_of_satellite(delta_lambda, lambda_0)
	m=float(lines[5])
	masses[k]=m
	print "Mass of molecules: ", m
	print "Velocity of satellite: ", velocity
	print "Temperature of gas:", compute_temperature_of_planet(float(lines[3]), lambda_0, m)
	print "Lambda_0:", 	lambda_0
	boundary_lambda=(0.1/600.0)*lambda_0
	lower_index_of_interest=np.argmax(spectrum[:, 0]>lambda_0-boundary_lambda)
	upper_index_of_interest=np.argmax(spectrum[:, 0]>lambda_0+boundary_lambda)
	interesting_lambdas=spectrum[lower_index_of_interest:upper_index_of_interest,0]
	interesting_flux=spectrum[lower_index_of_interest:upper_index_of_interest,1]
	guess=1+(float(lines[2])-1)*np.exp(-((interesting_lambdas-float(lines[4]))**2)/(2*(float(lines[3]))**2))
	plt.plot(interesting_lambdas, interesting_flux, label='_nolegend_')
	plt.plot(interesting_lambdas, guess)
	plt.title('Spectral line for ' + species[k])
	plt.xlabel('Wavelength [nm]')
	plt.ylabel('Flux [W/m^2]')
	plt.legend(['Best-fit line'])
	plt.show()
	happy=raw_input('Happy with this (y/n)? ')
	if happy.lower().startswith('y'):
		no_wavelengths[k]=True
	else:
		no_wavelengths[k]=False

	k+=1

existing_gasses=[]
previous_mass=masses[0]
prior_true=False
eps=1e-6
previously_existed=False
for i in range(len(masses)):
	if abs(previous_mass-masses[i])<eps:
		if no_wavelengths[i]==True:
			if previously_existed==True:
				pass
			else:
				existing_gasses.append(masses[i])
				previously_existed=True
		else:
			if previously_existed==True:
				if prior_true==True:
					del existing_gasses[-1]
			else:
				previously_existed=True
			prior_true=False
	else:
		previously_existed=False
		if no_wavelengths[i]==True:
			existing_gasses.append(masses[i])
		else:
			prior_true=False

print "Mean molar mass is: ", sum(existing_gasses)/float(len(existing_gasses))
			
				
				
		
