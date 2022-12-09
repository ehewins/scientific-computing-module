import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Input variables for Burrau's problem
G = 1 # Value of the gravitational constant
pos = np.array([[1,-2,1],[3,-1,-1]]) # initial positions of the three bodies
vel = np.array([[0,0,0],[0,0,0]]) # initial velocities of the three bodies
m = np.array([3,4,5]) # masses of the three bodies
num_dimensions = pos.shape[0]
num_bodies = pos.shape[1]
t_tot = 70 # total time to run the simulation for
dt = 0.001 # time step between outputs calculated by odeint
t_vals = np.arange(0,t_tot+dt,dt) # Array time values to calcuclate the state of the system at.
plot_llim = 60
plot_ulim = 70
in_range = (t_vals >= plot_llim) & (t_vals <= plot_ulim)
rtol = 1e-23 # Tolerance settings used for odeint. Default values are 1.49012e-8
atol = 1e-10

# # Input variables for the 2-body test case
# G = 6.6748015e-11
# pos = np.array([[-0.5*1.496e11,0.5*1.496e11],[0,0]])
# vel = np.array([[0,0],[-15*1e3,15*1e3]])
# m = np.array([1.989e30, 1.989e30])
# num_dimensions = pos.shape[0]
# num_bodies = pos.shape[1]
# t_tot = 5 * 31557600
# dt = 24*60*60
# t_vals = np.arange(0,t_tot+dt,dt)
# rtol = 1.49012e-8
# atol = 1.49012e-8

# # Input variables for two light bodies and a heavy body - Earth, Mars, Sun
# G = 6.6748015e-11
# pos = np.array([[-1.496e11,-1.524*1.496e11,0],[0,0,0]])
# vel = np.array([[0,0,0],[29.78*1e3,24.07*1e3,0]])
# m = np.array([5.9722e24, 6.4171e23, 1.989e30])
# num_dimensions = pos.shape[0]
# num_bodies = pos.shape[1]
# t_tot = 10 * 31557600
# dt = 24*60*60
# t_vals = np.arange(0,t_tot+dt,dt)
# rtol = 1.49012e-8
# atol = 1.49012e-8

def solve_func(X, t):
    '''Function for solving the differential equations of the problem, for n bodies in n dimensions.'''
    # NOTE: next thing to test is whether it's n dimensional yet.
    drdt = X[num_dimensions*num_bodies:]
    pos = X.reshape(2,num_dimensions,num_bodies)[0]
    dvdt = np.zeros((num_dimensions,num_bodies))
    for i in range(num_bodies):
        r_ij = pos - np.reshape(pos[:,i],(num_dimensions, 1))
        dists = np.sqrt(np.sum(r_ij**2, axis=0))
        dists[i] = 1 # prevents a divide by zero error
        dvdt[:,i] = G*np.sum(m*r_ij/dists**3, axis=1)
    return np.concatenate((drdt, dvdt), axis=None)

solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, rtol=rtol, atol=atol)

plt.figure(1)
plt.xlabel("$x$ coordinate")
plt.ylabel("$y$ coordinate")
for b in range(num_bodies):
    plt.plot(solutions[:,b][in_range], solutions[:,b+num_bodies][in_range], label='Body {:d}'.format(b+1))
plt.legend(loc='upper left')

total_energy = 0.5 * np.sum(np.tile(m,num_dimensions)*solutions[:,num_bodies*num_dimensions:]**2, axis=1) # First calculate total the kinetic energy, since this can be done on one line
for i in range(num_bodies): # Then add the gravitational potential energy by looping over each unique pair of bodies
    for j in range(i+1,num_bodies):
        total_energy += -G*m[i]*m[j]/np.sqrt(np.sum((solutions[:,j:num_bodies*num_dimensions:num_bodies]-solutions[:,i:num_bodies*num_dimensions:num_bodies])**2, axis=1))
Delta_E = (total_energy - total_energy[0]) / total_energy[0]
plt.figure(2)
plt.xlabel("$t$")
plt.ylabel("$\Delta E$")
plt.plot(t_vals, Delta_E)

plt.show()
