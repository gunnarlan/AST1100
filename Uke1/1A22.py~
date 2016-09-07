import numpy as np
import matplotlib.pyplot as plt

infile=np.loadtxt('Rainfall.dat')
##Method 1
monthly_mean_1=[]
for k in range(infile.shape[1]-2):
	 monthly_mean_1.append(np.mean(infile[:,k+1]))
print monthly_mean_1
##Method 2
monthly_mean_2=[]
nbins=20
for k in range(infile.shape[1]-2):
	largest=infile[:, k+1].max()
	smallest=infile[:, k+1].min()
	bins=np.linspace(smallest, largest, nbins+1)
	values=[]
	for l in range(len(bins)-1):
		relevant_values=np.logical_and(infile[:, k+1]>=bins[l], infile[:, k+1]<=bins[l+1])
		values.append(np.count_nonzero(relevant_values))
	probabilities=np.array(values)/(float(len(infile[:, k+1])))
	midpoints=[(a+b)/2.0 for a,b in zip(bins, bins[1::])]
	mean=np.sum(probabilities*midpoints)
	monthly_mean_2.append(mean)
print monthly_mean_2



