import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def grav_acc(i, r):
    '''Function which when passed an index i as well as an array of N position coordinates corresponding to all the bodies in the simulation and a second array of each body's mass, returns the net gravitational acceleration of the object at index i in the array'''
    G = 6.6748015e-11 # Value of the gravitational constant in Nm²kg⁻²
    r_ij = r - np.reshape(r[:,i],(num_dimensions, 1))
    r_ij = np.delete(r_ij, i, axis=1)
    M = np.delete(m, i)
    dists = np.sqrt(np.sum(r_ij**2, axis=0))
    return G*np.sum(M*r_ij/dists**3, axis=1)

def solve_func(X, t):
    X = np.reshape(X,(2,num_dimensions,num_bodies))
    pos = np.array(X[0])
    X[0] = X[1]
    for i in range(num_bodies):
        X[1,:,i] = grav_acc(i,pos)
    return X.flatten()

x1, x2 = -0.5 * 1.496e11, 0.5 * 1.496e11 # init. position components of bodies, in SI units
y1, y2 = 0, 0
vx1, vx2 = 0, 0 # init. velocity components of bodies, in SI units
vy1, vy2 = -15 * 1e3, 15 * 1e3
m1, m2 = 1.989e30, 1.989e30 # masses of the bodies, in SI units
pos = np.array([[x1,x2],[y1,y2]])
vel = np.array([[vx1,vx2],[vy1,vy2]])
m = np.array([m1,m2])
num_dimensions = np.shape(pos)[0]
num_bodies = np.shape(pos)[1]

t_tot = 5 * 31557600 # total time converted from years to seconds
dt = 24*60*60 # size of the timestep
t_vals = np.arange(0,t_tot+dt,dt) # Array of times
# print(np.concatenate((pos,vel),axis=None))
# print(solve_func(np.concatenate((pos,vel),axis=None), t_vals))
solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, full_output = 1)
x1_vals = solutions[:,0]
x2_vals = solutions[:,1]
y1_vals = solutions[:,2]
y2_vals = solutions[:,3]

fig = plt.figure(1,figsize=(14,5))
plt.subplots_adjust(wspace=0.25)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax1.plot(x1_vals/1.496e11,y1_vals/1.496e11,'r-',label='Body 1')
ax1.plot(x2_vals/1.496e11,y2_vals/1.496e11,'k-',label='Body 2')
ax1.legend(loc='upper left')
ax1.set_xlabel("$x$ [AU]")
ax1.set_ylabel("$y$ [AU]")
ax1.set_xlim(-.6,.6)
ax1.set_ylim(-.4,.4)
ax2.plot(t_vals/31557600, x1_vals/1.496e11, 'r-')
ax2.plot(t_vals/31557600, x2_vals/1.496e11, 'k-')
ax2.set_xlabel("$t$ [yr]")
ax2.set_ylabel("$x$ [AU]")

plt.show()
