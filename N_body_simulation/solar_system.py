import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.cm import tab20

def solve_func(X, t):
    '''
    Function used with scipy.integrate.odeint to solve the differential equations
    of the problem for n bodies in n dimensions.
    '''
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

def animate3d(frame):
    '''
    Function called by matplotlib.animation.FuncAnimation when animating a 3D system.
    '''
    body_index = 2*frame # we're only plotting a frame every 2nd datapoint
    body_size = 150
    tail_size = 1
    tail_length = 50 # Number of points to include in the tail which trails behind the main body
    if body_index <= tail_length:
        tail_end = 0
        tail_length = body_index
    else:
        tail_end = body_index - tail_length
    # Sets the positions of the scatter plot points
    bodies._offsets3d = np.array([
                                     solutions[tail_end:body_index+1,0:num_bodies].flatten(),
                                     solutions[tail_end:body_index+1,num_bodies:2*num_bodies].flatten(),
                                     solutions[tail_end:body_index+1,2*num_bodies:3*num_bodies].flatten()
                                 ])
    # Sets the sizes of the scatter plot points
    bodies.set_sizes(np.repeat(np.append(tail_size*np.ones(tail_length), [body_size]), num_bodies))
    # Sets the colours of the scatter plot points
    bodies.set_array(np.tile(np.arange(num_bodies), tail_length+1))
    return bodies

# Input variables for the sun and the 8 planets, with position and velocity data for 00:00 on 2000-01-01, obtained from NASA Horizons.
G = 6.6748015e-20 # Gravitational constant in km kg⁻¹(kms⁻¹)², reflecting the units of pos, vel and m
pos = np.array([[0,-2.105262107244070e7,-1.075055502719850e8,-2.521092855899356e7,2.079950549908331e8,5.989091645401344e8,9.587063371733198e8,2.158774699724352e9,2.514853282370434e9],[0,-6.640663812253430e7,-3.366520666522362e6,1.449279195838006e8,-3.143009561106971e6,4.391225866604841e8,9.825652104588115e8,-2.054825151185744e9,-3.738847414715512e9],[0,-3.492445946577720e6,6.159219789239045e6,-6.164165719002485e2,-5.178781160069674e6,-1.523251063025475e7,-5.522065631225652e7,-3.562361168065417e7,1.903959877100039e7]]) # Position coodinates of the sun + plantes in km
vel = np.array([[0,3.665298704187096e1,8.891597859686224e-1,-2.983983333677879e1,1.295003532851602e0,-7.901937516136118e0,-7.428885680409909e0,4.637648534301329e0,4.465682523947062e0],[0,-1.228983806940175e1,-3.515920774137907e1,-5.207633902410673e0,2.629442067068712e1,1.116317703172796e1,6.738814240733793e0,4.627193109110802e0,3.076493760667651e0],[0,-4.368173036243590e0,-5.318594228644749e-1,6.168441184239981e-5,5.190097267545717e-1,1.306732148714280e-1,1.776643606866641e-01,-4.285052612262108e-2,-1.657059897537549e-1]]) # Velocity components of the sun + plantes in km/s
m = np.array([1.988e30,3.302e23,4.8685e24,5.9722e24,6.4171e23,1.8982e27,5.6834e26,8.6813e25,1.02409e26]) # masses in kg
body_names = np.array(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"])

# # Simulating more bodies increases the time taken to perform the calculations, so if you do not wish to include all the planets in the simulation then uncomment this code block.
# # 0 is the Sun (don't omit this one), the planets are numbers 1 - 8 in their usual order.
# bodies_to_omit = np.arange(5,9) # only consider the inner solar system
# # bodies_to_omit = np.arange(1,5) # only consider the outer solar system
# # bodies_to_omit = np.array([6,7,8]) # custom selection
# pos = np.delete(pos, bodies_to_omit, axis=1)
# vel = np.delete(vel, bodies_to_omit, axis=1)
# m = np.delete(m, bodies_to_omit)
# body_names = np.delete(body_names, bodies_to_omit)

# This section is used to choose an appropriate length of time to run the simulation for.
# For rogue body close encounters which simulate only inner solar system, longer timescalse may be required depending on the initial position/velocity of the body
expected_orbital_periods = np.array([0, 88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800])*24*60*60 # Expected orbital periods of the bodies, converted from days to seconds.
try: # This will fail if bodies_to_omit is commented out, in which case the simulation assumes the required total time is the orbital period of Neptune.
    t_tot = np.delete(expected_orbital_periods, bodies_to_omit)[-1]
except:
    t_tot = expected_orbital_periods[-1]

# t_tot = t_tot * 20 # Use this line to adjust if the time given by the above block is insufficient

dt = 24*60*60 # timesteps at which odeint should calculate positions and velocities
t_vals = np.arange(0,t_tot+dt,dt)
rtol = 1e-4 # integration tolerances used by odeint, increased for speed. Default: 1.49012e-8
atol = 1e-4

# # The next seven lines may be uncommented if you want to include an additional body in the solar system simulation, such as a rogue star.
# newbody_pos = np.array([10*1.5e8, -30*1.5e8, 0])
# # newbody_vel = np.array([0, 2*10*1.5e8/t_tot, 0])
# newbody_vel = np.array([0, 10, 0]) # rogue stars tend to have velocities on the order of 10km/s
# newbody_mass = 0.5*1.988e30
# pos = np.append(pos, np.array([newbody_pos]).T, axis=1)
# vel = np.append(vel, np.array([newbody_vel]).T, axis=1)
# m = np.append(m, newbody_mass)
# body_names = np.append(body_names, "Rogue body")

num_dimensions = pos.shape[0]
num_bodies = pos.shape[1]

fix_first_body = False # Set this variable depending on whether you want body 0 (the Sun) to be affected by the gravity of the planets. Warning: violates energy conservation (but generates prettier graphs)
print("Beginning calculations...")
start_time = time.time()
solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, rtol=rtol, atol=atol)
end_time = time.time()
print("Solving the differential equations took",end_time-start_time,"seconds")

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
for b in range(num_bodies):
    ax.plot(solutions[:,b], solutions[:,b+num_bodies], solutions[:,b+2*num_bodies], label=body_names[b])
ax.legend(loc='upper left', fontsize=9, ncol=2 if num_bodies>5 else 1)
ax.set_xlabel("$x$ position [km]")
ax.set_ylabel("$y$ position [km]")
ax.set_zlabel("$z$ position [km]")

# Calculating the relative change in the total energy of the system:
# First calculate total the kinetic energy, since this can be done on one line
total_energy = 0.5 * np.sum(np.tile(m,num_dimensions)*solutions[:,num_bodies*num_dimensions:]**2, axis=1)
# Then add the gravitational potential energy by looping over each unique pair of bodies
for i in range(num_bodies):
    for j in range(i+1,num_bodies):
        total_energy += -G*m[i]*m[j]/np.sqrt(np.sum((solutions[:,j:num_bodies*num_dimensions:num_bodies]-solutions[:,i:num_bodies*num_dimensions:num_bodies])**2, axis=1))
Delta_E = (total_energy - total_energy[0]) / total_energy[0]
plt.figure(2)
plt.xlabel("Time, $t$")
plt.ylabel("Relative change in total energy, $\Delta E$")
plt.plot(t_vals, Delta_E)

print("Would you like to generage an animation of the simulation you just ran? [y/n]")
animate_yn = input("Your choice: ").lower()
if animate_yn in ("y", "yes"):
    anifig = plt.figure(3)
    aniax = anifig.add_subplot(projection='3d')
    bodies = aniax.scatter([], [], [], c=[], marker='.', cmap='tab20')
    anim = FuncAnimation(anifig, animate3d, frames=t_vals.size//2, interval=70/t_tot)
    aniax.set_xlim(1.05*np.min(solutions[:,0:num_bodies]), 1.05*np.max(solutions[:,0:num_bodies]))
    aniax.set_ylim(1.05*np.min(solutions[:,num_bodies:2*num_bodies]), 1.05*np.max(solutions[:,num_bodies:2*num_bodies]))
    aniax.set_zlim(1.05*np.min(solutions[:,2*num_bodies:3*num_bodies]), 1.05*np.max(solutions[:,2*num_bodies:3*num_bodies]))
    aniax.set_xlabel("$x$ coordinate")
    aniax.set_ylabel("$y$ coordinate")
    aniax.set_zlabel("$z$ coordinate")
else:
    print("You have chosen not to generate an animation.")

plt.show()
