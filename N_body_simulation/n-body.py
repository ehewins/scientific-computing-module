import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def grav_acc(i, r):
    '''Function which when passed an index i as well as an array of N position coordinates corresponding to all the bodies in the simulation, returns the net gravitational acceleration of the object at index i in the array'''
    G = 6.6748015e-11 # Value of the gravitational constant in Nm²kg⁻²
    r_ij = r - np.reshape(r[:,i],(num_dimensions, 1))
    r_ij = np.delete(r_ij, i, axis=1)
    M = np.delete(m, i)
    dists = np.sqrt(np.sum(r_ij**2, axis=0))
    return G*np.sum(M*r_ij/dists**3, axis=1)

plt.show()
