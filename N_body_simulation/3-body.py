"""
This script is another test, which aims to simulate three bodies in 2D rather than just two. The goal is to demonstrate the solution to Burrau's problem shown in the paper by Szebehely and Peters, DOI: 10.1086/110355
"""

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
t_tot = 10 # total time to run the simulation for
dt = 0.001 # time step between outputs calculated by odeint
t_vals = np.arange(0,t_tot+dt,dt) # Array time values to calcuclate the state of the system at.
rtol = 1.49012e-10 # Tolerance settings used for odeint. Default values are 1.49012e-8
atol = 1.49012e-10

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

# def solve_func(X, t):
#     '''Function for solving the differential equations of the problem, for 3 bodies in 2 dimensions.'''
#     x1, x2, x3, y1, y2, y3, vx1, vx2, vx3, vy1, vy2, vy3 = X
#     m1, m2, m3 = m
#     dx1dt, dx2dt, dx3dt = vx1, vx2, vx3
#     dy1dt, dy2dt, dy3dt = vy1, vy2, vy3
#     dvx1dt = G*m2*(x2-x1)/np.sqrt((x2-x1)**2+(y2-y1)**2)**3 + G*m3*(x3-x1)/np.sqrt((x3-x1)**2+(y3-y1)**2)**3
#     dvy1dt = G*m2*(y2-y1)/np.sqrt((x2-x1)**2+(y2-y1)**2)**3 + G*m3*(y3-y1)/np.sqrt((x3-x1)**2+(y3-y1)**2)**3
#     dvx2dt = G*m1*(x1-x2)/np.sqrt((x1-x2)**2+(y1-y2)**2)**3 + G*m3*(x3-x2)/np.sqrt((x3-x2)**2+(y3-y2)**2)**3
#     dvy2dt = G*m1*(y1-y2)/np.sqrt((x1-x2)**2+(y1-y2)**2)**3 + G*m3*(y3-y2)/np.sqrt((x3-x2)**2+(y3-y2)**2)**3
#     dvx3dt = G*m1*(x1-x3)/np.sqrt((x1-x3)**2+(y1-y3)**2)**3 + G*m2*(x2-x3)/np.sqrt((x2-x3)**2+(y2-y3)**2)**3
#     dvy3dt = G*m1*(y1-y3)/np.sqrt((x1-x3)**2+(y1-y3)**2)**3 + G*m2*(y2-y3)/np.sqrt((x2-x3)**2+(y2-y3)**2)**3
#     return [dx1dt, dx2dt, dx3dt, dy1dt, dy2dt, dy3dt, dvx1dt, dvx2dt, dvx3dt, dvy1dt, dvy2dt, dvy3dt]
    
def solve_func(X, t):
    '''Function for solving the differential equations of the problem, for n bodies in 2 dimensions.'''
    drdt = X[num_dimensions*num_bodies:]
    pos = X.reshape(2,2,num_bodies)[0]
    dvxdt = np.zeros(num_bodies)
    dvydt = np.zeros(num_bodies)
    for i in range(num_bodies):
        for j in range(num_bodies):
            if j != i:
                distance = np.sqrt((pos[0,j]-pos[0,i])**2 + (pos[1,j]-pos[1,i])**2)
                dvxdt[i] += G*m[j]*(pos[0,j]-pos[0,i])/distance**3
                dvydt[i] += G*m[j]*(pos[1,j]-pos[1,i])/distance**3
    return np.concatenate((drdt, dvxdt, dvydt))

solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, rtol=rtol, atol=atol)

plt.figure(1)
plt.xlabel("$x$ coordinate")
plt.ylabel("$y$ coordinate")
for b in range(num_bodies):
    plt.plot(solutions[:,b], solutions[:,b+num_bodies])

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
