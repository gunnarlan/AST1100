import numpy as np
from AST1100SolarSystem import AST1100SolarSystem

class compute_angle:
	def __init__(self, seed):
		self.mysolarsystem=AST1100SolarSystem(seed,hasMoons=False)
		self.target_planet=6
		self.angular_vel=2.0*np.pi/(60*60.0*24*self.mysolarsystem.period[self.target_planet])
		self.period_of_satellite=6899.30546952

	def find_angle_of_satellite(self, position, picture_times):
		r=np.linalg.norm(position, axis=0)
		position=-position
		phi=np.arctan2(position[1,:], position[0,:])
		theta=np.arccos(position[2,:]/(1.0*r))
		return theta, phi



infile=open('Pictures.txt', 'r') #File containing the ouput from the orientation of the satellite
lines=infile.readlines()
no_of_pictures=21
t0=171000
period_of_satellite=6899.30546952
picture_times=np.linspace(t0, period_of_satellite+t0, no_of_pictures)
positions=np.zeros(shape=(3, no_of_pictures))
j=0
for k in range(len(lines)):
	if lines[k].startswith('Pos'):
		values=[]
		coordinates=lines[k].split('[')
		coordinates= coordinates[1].split(']')[0]
		coordinates=coordinates.split(' ')
		for i in range(len(coordinates)):
			try:
				float(coordinates[i])
				values.append(float(coordinates[i]))
			except ValueError:
				pass
		positions[:, j]=values
		j+=1

comp=compute_angle(6392)
theta, phi=comp.find_angle_of_satellite(positions, picture_times)
for k in range(len(theta)):
	print "picture "+str(picture_times[k])+' '+str(theta[k]) + ' ' +str(phi[k]) + ' 0 0 1'

