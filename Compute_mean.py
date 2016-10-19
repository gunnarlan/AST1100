import numpy as np
import itertools
possible_masses=[32,18.02,44.01,16.04,28.01,44.02]

for L in range(0, len(possible_masses)+1):
	for subset in itertools.combinations(possible_masses, L):
		print subset
		print np.mean(subset)


