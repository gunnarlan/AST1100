import numpy as np
import matplotlib.pyplot as plt

infile=np.loadtxt('Rainfall.dat')
january=infile[:,1]
nbins=20
largest=january.max()
smallest=january.min()
bins=np.linspace(smallest, largest, nbins+1)
values=[]
for k in range(len(bins)-1):
	relevant_values=np.logical_and(january>=bins[k], january<=bins[k+1])
	values.append(np.count_nonzero(relevant_values))

values=np.array(values)/(float(len(january)))
midpoints=[(a+b)/2.0 for a,b in zip(bins, bins[1::])]
print len(midpoints)
plt.plot(midpoints, values)
plt.title('Histogram for January')
plt.xlabel('Amount of rainfall')
plt.ylabel('Probability')
plt.show()

