import numpy as np
import matplotlib.pyplot as plt

def Gauss(x,sigma, mu):
	prefac=1/(sigma*np.sqrt(2*np.pi))
	exponent=-((x-mu)**2)/(2.0*sigma**2)
	return prefac*np.exp(exponent)

infile=np.loadtxt('Rainfall.dat')
yearly_means=[]
for l in range(infile.shape[0]):
	current_means=[]
	for k in range(infile.shape[1]-2):
		current_means.append(infile[l][k+1])
	yearly_means.append(np.mean(current_means))
	
yearly_means=np.array(yearly_means)

mean_of_means=np.mean(yearly_means)
std_of_means=np.std(yearly_means)


nbins=20
largest=yearly_means.max()
smallest=yearly_means.min()
x=np.linspace(smallest, largest, 10000)
bins=np.linspace(smallest, largest, nbins+1)
values=[]
for k in range(len(bins)-1):
	relevant_values=np.logical_and(yearly_means>=bins[k], yearly_means<=bins[k+1])
	values.append(np.count_nonzero(relevant_values))

probabilities=np.array(values)/(float(len(yearly_means)))
midpoints=[(a+b)/2.0 for a,b in zip(bins, bins[1::])]

figure1=plt.figure()
plt.plot(x, Gauss(x, std_of_means, mean_of_means))
plt.plot(midpoints, probabilities)
plt.title('Histogram for yearly means')
plt.xlabel('Amount of rainfall')
plt.ylabel('Probability')
plt.show()
