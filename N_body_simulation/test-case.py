"""
This script simulates a simple test where there are only two bodies of equal mass orbiting around each other in two dimensions. It solves 8 coupled differential equations for x & y components of position & velocity, that's four equations for each of the bodies. The code is not generalisable in its current form, but this it serves to establish the basic principles involved in the project.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def solve_func(X,t,m1,m2):
    G = 6.6748015e-11 # G is in SI units, so all input parameters should also be in SI units.
    # G = 8.868e24 # G in astronomical units per solar mass kilometers squared per second squared
    x1,x2,y1,y2,vx1,vx2,vy1,vy2 = X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7]
    distance = np.sqrt((x2-x1)**2+(y2-y1)**2)
    # The 8 differential equations which needed solving
    dx1dt = vx1
    dy1dt = vy1
    dvx1dt = G*m2*(x2-x1)/distance**3
    dvy1dt = G*m2*(y2-y1)/distance**3
    dx2dt = vx2
    dy2dt = vy2
    dvx2dt = G*m1*(x1-x2)/distance**3
    dvy2dt = G*m1*(y1-y2)/distance**3
    return [dx1dt, dx2dt, dy1dt, dy2dt, dvx1dt, dvx2dt, dvy1dt, dvy2dt]

x1, y1 = -0.5 * 1.496e11, 0 # init. position components of bodies, in SI units
x2, y2 = 0.5 * 1.496e11, 0
vx1, vy1 = 0, -15 * 1e3 # init. velocity components of bodies, in SI units
vx2, vy2 = 0, 15 * 1e3
m1, m2 = 1.989e30, 1.989e30 # masses of the bodies, in SI units
t_tot = 5 * 31557600 # total time converted from years to seconds
dt = 24*60*60 # size of the timestep

t_vals = np.arange(0,t_tot+dt,dt) # Array of times
# print(solve_func((x1,x2,y1,y2,vx1,vx2,vy1,vy2),t_vals,m1,m2))
solutions = odeint(solve_func, (x1,x2,y1,y2,vx1,vx2,vy1,vy2), t_vals, args=(m1,m2)) # Call to odeint to slove the equations
x1_vals = solutions[:,0] # Extracting the solutions for the positions of the bodies across time
x2_vals = solutions[:,1]
y1_vals = solutions[:,2]
y2_vals = solutions[:,3]

# Displaying a figure to prove the test case was successful.
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

# Calculating the relative change in total energy across time, and plotting it to ensure energy is conserved.
G = 6.6748015e-11
total_gravitational = -G*m1*m2/np.sqrt((solutions[:,1]-solutions[:,0])**2 + (solutions[:,3]-solutions[:,2])**2)
total_kinetic = 0.5*m1*(solutions[:,4]**2 + solutions[:,6]**2) + 0.5*m2*(solutions[:,5]**2 + solutions[:,7]**2)
total_energy = total_gravitational + total_kinetic
Delta_E = (total_energy - total_energy[0]) / total_energy[0]
plt.figure(2)
plt.plot(t_vals/31557600,Delta_E)
plt.xlabel("$t$ [yr]")
plt.ylabel("$\Delta E$")
plt.title("Relative change in the total energy with time.\nIf $|\Delta E|<10^{-5}$, the code conserves energy.")

plt.show()
