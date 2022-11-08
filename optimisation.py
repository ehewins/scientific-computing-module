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
# deliberately slow code to find nearest neighbours within periodic unit cube
#
#  You may only change the code between here and the line "end_time=time.time()")
#
s = np.zeros(shape=(N,N))
matchedIndices = np.zeros(N)
for a in range(N):
    for b in range(N):
        distance_vectors = np.abs(pos[:,a] - pos[:,b])
        s[a,b] = np.sum(np.minimum(distance_vectors, 1-distance_vectors)**2)

s = np.sqrt(s)

for a in range(N):
    mindist = 1e10
    matchedIndices[a] = 0
    for c in range(N):
        if a != c:
            mindist = np.minimum(s[a,c], mindist)
            if mindist == s[a,c]:
                matchedIndices[a] = c

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
