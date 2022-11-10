# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 12:00:35 2016

@author: ppzfrp, modified by ppyeh5
"""
# set up random 3-d positions
#
import numpy as np
import time
import errno
N = 100
seed = 1234
np.random.seed(seed)
pos = np.random.random((3,N))
start_time = time.time()
# deliberately slow code to find nearest neighbours within periodic unit cube - It's not slow now that I've finished with it!
#
#  You may only change the code between here and the line "end_time=time.time()")
#
matchedIndices = np.zeros(N)
for a in range(N):
    dists_array = np.abs(pos - np.reshape(pos[:,a],(3,1))) # Creates a (3,N) array of the distances between the coordinate given by column a, and all the other columns.
    min_dists = np.minimum(dists_array, 1-dists_array) # Finds the actual sortest distances, given the periodic boundaries of the unit cube.
    s = np.sum(min_dists**2, axis=0) # Produces an array of the square magnitudes of the displacement beteen a and all other points by summing the squares of each column.
    s[a] = 1 # Sets the distance from point a to itself to one, so that it won't get counted as its own nearest neighbor.
    matchedIndices[a] = np.argmin(s) # Saves the index of point a's nearest neighbor.

end_time = time.time()
print('Elapsed time = ', repr(end_time - start_time))

# generate filename from N and seed
filename = 'pyneigh' + str(N) + '_' + str(seed)
# if a file with this name already exists read in the nearest neighbour
# list found before and compare this to the current list of neighbours,
# else save the neighbour list to disk for future reference
try:
    fid = open(filename,'rb')
    matchedIndicesOld = np.loadtxt(fid)
    fid.close()
    if (matchedIndicesOld == matchedIndices).all():
        print('Checked match')
    else:
        print('Failed match')
except OSError as e:
    if e.errno == errno.ENOENT:
        print('Saving neighbour list to disk')
        fid = open(filename,'wb')
        np.savetxt(fid, matchedIndices, fmt="%8i")
        fid.close()
    else:
        raise
