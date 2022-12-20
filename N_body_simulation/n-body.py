import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.cm import tab20

fix_first_body = False

# # Input variables for Burrau's problem
# G = 1 # Value of the gravitational constant
# pos = np.array([[1,-2,1],[3,-1,-1]]) # initial positions of the three bodies
# vel = np.array([[0,0,0],[0,0,0]]) # initial velocities of the three bodies
# m = np.array([3,4,5]) # masses of the three bodies
# num_dimensions = pos.shape[0]
# num_bodies = pos.shape[1]
# t_tot = 60 # total time to run the simulation for
# dt = 0.001 # time step between outputs calculated by odeint
# t_vals = np.arange(0,t_tot+dt,dt) # Array time values to calcuclate the state of the system at.
# plot_llim = 50
# plot_ulim = 60
# # rtol = 1.0001e-24 # Tolerance settings used for odeint. Default values are 1.49012e-8
# rtol = 1e-26 # Tolerance settings used for odeint. Default values are 1.49012e-8
# # atol = 8.114e-11
# atol = 9.549925860215698e-11 # Best energy conserver (min. weighted average method)

# Input variables for simulation of the solar system, from NASA Horizons at 00:00 on 2000-01-01
G = 6.6748015e-20 # Gravitational constant in km kg⁻¹(kms⁻¹)²
pos = np.array([[0,-2.105262107244070e7,-1.075055502719850e8,-2.521092855899356e7,2.079950549908331e8,5.989091645401344e8,9.587063371733198e8,2.158774699724352e9,2.514853282370434e9],[0,-6.640663812253430e7,-3.366520666522362e6,1.449279195838006e8,-3.143009561106971e6,4.391225866604841e8,9.825652104588115e8,-2.054825151185744e9,-3.738847414715512e9],[0,-3.492445946577720e6,6.159219789239045e6,-6.164165719002485e2,-5.178781160069674e6,-1.523251063025475e7,-5.522065631225652e7,-3.562361168065417e7,1.903959877100039e7]]) # Position coodinates of the sun + plantes in km
vel = np.array([[0,3.665298704187096e1,8.891597859686224e-1,-2.983983333677879e1,1.295003532851602e0,-7.901937516136118e0,-7.428885680409909e0,4.637648534301329e0,4.465682523947062e0],[0,-1.228983806940175e1,-3.515920774137907e1,-5.207633902410673e0,2.629442067068712e1,1.116317703172796e1,6.738814240733793e0,4.627193109110802e0,3.076493760667651e0],[0,-4.368173036243590e0,-5.318594228644749e-1,6.168441184239981e-5,5.190097267545717e-1,1.306732148714280e-1,1.776643606866641e-01,-4.285052612262108e-2,-1.657059897537549e-1]]) # Velocity components of the sun + plantes in km/s
m = np.array([1.988e30,3.302e23,4.8685e24,5.9722e24,6.4171e23,1.8982e27,5.6834e26,8.6813e25,1.02409e26]) # masses in kg
num_dimensions = pos.shape[0]
num_bodies = pos.shape[1]
# t_tot = 164.8 * 31557600 # chosen to be the orbital period of neptune in seconds
t_tot = 2* 164.8 * 31557600 # chosen to be the orbital period of neptune in seconds
dt = 24*60*60 * 7
t_vals = np.arange(0,t_tot+dt,dt)
plot_llim = 0
plot_ulim = t_tot
rtol = 1.49012e-8
atol = 1.49012e-8
# fix_first_body = True
# The next seven lines may be uncommented if you want to include an additional body in the solar system simulation, such as a rouge star.
newbody_pos = np.array([80*1.5e8, -80*1.5e8, 0])
newbody_vel = np.array([0, 2*80*1.5e8/t_tot, 0])
newbody_mass = 0.5*1.988e30
pos = np.append(pos, np.array([newbody_pos]).T, axis=1)
vel = np.append(vel, np.array([newbody_vel]).T, axis=1)
m = np.append(m, newbody_mass)
num_bodies += 1

# # Input variables for a 3 dimensional four body system
# G = 1
# m = np.array([1000,1,1,1]) # The central body should have a significantly greater mass than the other bodies in the system
# pos = np.array([[0,10,0,0],[0,0,20,0],[0,0,0,30]])
# # For this case, the velocity values are calculated based on the positions so as to give stable orbits
# vel = np.array([[0,0,0,np.sqrt(G*m[0]/pos[2,3])],[0,np.sqrt(G*m[0]/pos[0,1]),0,0],[0,0,np.sqrt(G*m[0]/pos[1,2]),0]])
# num_dimensions = pos.shape[0]
# num_bodies = pos.shape[1]
# t_tot = 50
# dt = 0.001
# t_vals = np.arange(0,t_tot+dt,dt)
# plot_llim = 0
# plot_ulim = t_tot
# rtol = 1.49012e-8
# atol = 1.49012e-8

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
# plot_llim = 0
# plot_ulim = t_tot
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
# plot_llim = 0
# plot_ulim = t_tot
# rtol = 1.49012e-8
# atol = 1.49012e-8

def solve_func(X, t):
    '''Function for solving the differential equations of the problem, for n bodies in n dimensions.'''
    drdt = X[num_dimensions*num_bodies:]
    pos = X.reshape(2,num_dimensions,num_bodies)[0]
    dvdt = np.zeros((num_dimensions,num_bodies))
    for i in range(num_bodies):
        r_ij = pos - np.reshape(pos[:,i],(num_dimensions, 1))
        dists = np.sqrt(np.sum(r_ij**2, axis=0))
        dists[i] = 1 # prevents a divide by zero error. r_ij[i] = 0 anyway, so it doesn't affect the next line
        dvdt[:,i] = G*np.sum(m*r_ij/dists**3, axis=1)
    # If the first body is at the centre of the orbital system, then these two lines will ensure it remains at the centre of the coordinate system.
    if fix_first_body == True:
        drdt[::num_bodies] = 0
        dvdt[:,0] = 0
    return np.concatenate((drdt, dvdt), axis=None)

start_time = time.time()
solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, rtol=rtol, atol=atol)
end_time = time.time()
print("Solving the differential equations took",end_time-start_time,"seconds")

fig = plt.figure(1)
in_range = (t_vals >= plot_llim) & (t_vals <= plot_ulim)
if num_dimensions == 2:
    ax = fig.add_subplot()
    for b in range(num_bodies):
        ax.plot(solutions[:,b][in_range], solutions[:,b+num_bodies][in_range], label='Body {:d}'.format(b+1))
elif num_dimensions == 3:
    ax = fig.add_subplot(projection='3d')
    for b in range(num_bodies):
        ax.plot(solutions[:,b][in_range], solutions[:,b+num_bodies][in_range], solutions[:,b+2*num_bodies][in_range], label='Body {:d}'.format(b+1))
    ax.set_zlabel("$z$ coordinate")
else:
    print("Improper number of dimensions, must be 2D or 3D for data to be plotted.")
    exit()

ax.legend(loc='upper left')
ax.set_xlabel("$x$ coordinate")
ax.set_ylabel("$y$ coordinate")

total_energy = 0.5 * np.sum(np.tile(m,num_dimensions)*solutions[:,num_bodies*num_dimensions:]**2, axis=1) # First calculate total the kinetic energy, since this can be done on one line
for i in range(num_bodies): # Then add the gravitational potential energy by looping over each unique pair of bodies
    for j in range(i+1,num_bodies):
        total_energy += -G*m[i]*m[j]/np.sqrt(np.sum((solutions[:,j:num_bodies*num_dimensions:num_bodies]-solutions[:,i:num_bodies*num_dimensions:num_bodies])**2, axis=1))
Delta_E = (total_energy - total_energy[0]) / total_energy[0]
plt.figure(2)
plt.xlabel("$t$")
plt.ylabel("$\Delta E$")
plt.plot(t_vals, Delta_E)
print("Calculated with a relative error range of",np.abs(np.max(Delta_E)-np.min(Delta_E)))

def animate2d(frame):
    body_index = 10*frame
    tail_length = 250
    body_size = 150
    tail_size = 0.2 # number of datapoints to include in the "tail" which trails behind the main body.
    if body_index <= tail_length:
        tail_end = 0
        tail_length = body_index
    else:
        tail_end = body_index - tail_length
    bodies.set_offsets(solutions[tail_end:body_index+1,:2*num_bodies].reshape((num_bodies*(tail_length+1),2), order='F'))
    bodies.set_sizes(np.append(tail_size*np.ones(tail_length), [body_size]))
    bodies.set_array(np.repeat(np.arange(num_bodies), tail_length+1))

def animate3d(frame):
    body_index = 10*frame
    tail_length = 800 # number of datapoints to include in the "tail" which trails behind the main body.
    body_size = 150
    tail_size = 0.2
    if body_index <= tail_length:
        tail_end = 0
        tail_length = body_index
    else:
        tail_end = body_index - tail_length
    bodies._offsets3d = np.array([
                                     solutions[tail_end:body_index+1,0:num_bodies].flatten(),
                                     solutions[tail_end:body_index+1,num_bodies:2*num_bodies].flatten(),
                                     solutions[tail_end:body_index+1,2*num_bodies:3*num_bodies].flatten()
                                 ])
    bodies.set_sizes(np.repeat(np.append(tail_size*np.ones(tail_length), [body_size]), num_bodies))
    bodies.set_array(np.tile(np.arange(num_bodies), tail_length+1))

# anifig = plt.figure(3)
# if num_dimensions == 2:
#     aniax = anifig.add_subplot()
#     bodies = aniax.scatter([], [], c=[], marker='.', cmap='tab20')
#     anim = FuncAnimation(anifig, animate2d, frames=t_vals.size//10, interval=20)
# elif num_dimensions == 3:
#     aniax = anifig.add_subplot(projection='3d')
#     bodies = aniax.scatter([], [], [], c=[], marker='.', cmap='tab20')
#     anim = FuncAnimation(anifig, animate3d, frames=t_vals.size//10, interval=1e-5)
#     aniax.set_zlim(1.05*np.min(solutions[:,2*num_bodies:3*num_bodies]), 1.05*np.max(solutions[:,2*num_bodies:3*num_bodies]))
#     aniax.set_zlabel("$z$ coordinate")
# else:
#     print("Improper number of dimensions, must be 2D or 3D for data to be plotted.")
#     exit()
# aniax.set_xlim(1.05*np.min(solutions[:,0:num_bodies]), 1.05*np.max(solutions[:,0:num_bodies]))
# aniax.set_ylim(1.05*np.min(solutions[:,num_bodies:2*num_bodies]), 1.05*np.max(solutions[:,num_bodies:2*num_bodies]))
# aniax.set_xlabel("$x$ coordinate")
# aniax.set_ylabel("$y$ coordinate")

plt.show()
