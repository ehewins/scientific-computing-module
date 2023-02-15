"""
The goal of this script is to find the integration tolerances, rtol and atol, which lead to the most accurate solution to Burrau's problem, replicating the results shown in the paper by Szebehely and Peters, DOI: 10.1086/110355
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation

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
# rtol = 1.0001e-24 # Tolerance settings used for odeint. Default values are 1.49012e-8
rtol = 1e-26 # Tolerance settings used for odeint. Default values are 1.49012e-8
# Because the values of atol will be so small, it makes sence to space the different values logarithmically rather than linearly. Only varying atol because it seems to have a more significant effect on the simulation than rtol.
log_atol = np.arange(-10.2,-9.8,0.01)
# log_atol = np.arange(-10.2,-10,0.002)
atol_vals = 10**log_atol
atol = 8.114e-11 # This is an initial value of atol, based on the best values found via trial and error. From the energy graph this produces, a model for good subsiquent values will be generated.

def solve_func(X, t):
    '''Function for solving the differential equations of the problem, for n bodies in n dimensions.'''
    # NOTE: next thing to test is whether it's n dimensional yet.
    drdt = X[num_dimensions*num_bodies:]
    pos = X.reshape(2,num_dimensions,num_bodies)[0]
    dvdt = np.zeros((num_dimensions,num_bodies))
    for i in range(num_bodies):
        r_ij = pos - np.reshape(pos[:,i],(num_dimensions, 1))
        dists = np.sqrt(np.sum(r_ij**2, axis=0))
        dists[i] = 1 # prevents a divide by zero error. r_ij[i] = 0 anyway, so it doesn't affect the next line
        dvdt[:,i] = G*np.sum(m*r_ij/dists**3, axis=1)
    return np.concatenate((drdt, dvdt), axis=None)

print("Calculating initial solution")
solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, rtol=rtol, atol=atol)
total_energy = 0.5 * np.sum(np.tile(m,num_dimensions)*solutions[:,num_bodies*num_dimensions:]**2, axis=1) # First calculate total the kinetic energy, since this can be done on one line
for i in range(num_bodies): # Then add the gravitational potential energy by looping over each unique pair of bodies
    for j in range(i+1,num_bodies):
        total_energy += -G*m[i]*m[j]/np.sqrt(np.sum((solutions[:,j:num_bodies*num_dimensions:num_bodies]-solutions[:,i:num_bodies*num_dimensions:num_bodies])**2, axis=1))
Delta_E = total_energy - total_energy[0]
# Finding the convolution smooths out the spikes, as they will be in slightly different locations for different plots. As such, they should not affect the weighting
smoothing_param = 150
smoothed = np.convolve(Delta_E, np.ones(smoothing_param)/smoothing_param, mode='same')
smoothed[np.argmax(smoothed):] = np.max(smoothed) # Without this line, the graph of the smoothed data suddenly drops towards the end, which we don't want

# Calculates an array of the average change in total energy during the simulation, weighted so as to emphasise energy conservation where it is most important
print("Beginning search for optimal atol values")
start_time = time.time()
err_avg = np.array([])
for atol in atol_vals:
    solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, rtol=rtol, atol=atol)
    total_energy = 0.5 * np.sum(np.tile(m,num_dimensions)*solutions[:,num_bodies*num_dimensions:]**2, axis=1) # First calculate total the kinetic energy, since this can be done on one line
    for i in range(num_bodies): # Then add the gravitational potential energy by looping over each unique pair of bodies
        for j in range(i+1,num_bodies):
            total_energy += -G*m[i]*m[j]/np.sqrt(np.sum((solutions[:,j:num_bodies*num_dimensions:num_bodies]-solutions[:,i:num_bodies*num_dimensions:num_bodies])**2, axis=1))
    Delta_E = total_energy - total_energy[0]
    err_avg = np.append(err_avg, np.abs(np.average(Delta_E, weights=smoothed)))
end_time = time.time()
print("Solving the equations for all atol values took",end_time-start_time,"seconds")

plt.figure(1)
ax1 = plt.axes()
ax1.set_xlabel("Value of atol")
ax1.set_ylabel("Weighted average of relative change in energy")
ax1.plot(atol_vals[err_avg<1e-5], err_avg[err_avg<1e-5], '.-') # Some of the energy values are >> 1e-5, so there's no point plotting them since they're obviously wrong.
plt.figure(2)
ax2 = plt.axes()
ax2.set_xlabel("log10 of atol")
ax2.set_ylabel("Weighted average of relative change in energy")
ax2.plot(np.log10(atol_vals[err_avg<1e-5]), err_avg[err_avg<1e-5], '.-')
# Print the best energy conserving solutions, in order of minimum energy change first
for i in range(atol_vals.size):
    print("Err_avg of",np.sort(err_avg)[i],"due to atol =",atol_vals[np.argsort(err_avg)][i])
print("\nHighest accuracy found using atol =",atol_vals[np.argmin(err_avg)])

plt.show()
