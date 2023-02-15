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
    return np.concatenate((drdt, dvdt), axis=None)

def animate2d(frame):
    '''
    Function called by matplotlib.animation.FuncAnimation when animating a 2D system.
    '''
    body_index = 10*frame # we're only plotting a frame every 10th datapoint
    body_size = 150
    tail_size = 0.2
    tail_length = animation_tail_length # Number of points to include in the tail which trails behind the main body
    if body_index <= tail_length:
        tail_end = 0
        tail_length = body_index
    else:
        tail_end = body_index - tail_length
    # Sets the positions of the scatter plot points
    bodies.set_offsets(solutions[tail_end:body_index+1,:2*num_bodies].reshape((num_bodies*(tail_length+1),2), order='F'))
    # Sets the sizes of the scatter plot points
    bodies.set_sizes(np.append(tail_size*np.ones(tail_length), [body_size]))
    # Sets the colours of the scatter plot points
    bodies.set_array(np.repeat(np.arange(num_bodies), tail_length+1))
    return bodies

def animate3d(frame):
    '''
    Function called by matplotlib.animation.FuncAnimation when animating a 3D system.
    '''
    body_index = 10*frame # we're only plotting a frame every 10th datapoint
    body_size = 150
    tail_size = 0.2
    tail_length = animation_tail_length # Number of points to include in the tail which trails behind the main body
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

# User input section for selecting which simulation to run
print("Please choose one of the following simulations to run, and enter its number.")
print("1 - Test case: 2 bodies in two dimensions, orbiting one another.")
print("2 - 3 body system: Two light bodies orbiting a much more massive one in 2 dimensions")
print("3 - 3D simulation showcase: three light bodies orbiting a heavy body in 3 dimensions")
print("4 - Burrau's problem in 2D")
print("5 - Burrau's problem in 3D")
choice = int(input("Your choice: "))

# Settings which are the same for the majority of the simulations
rtol = 1.49012e-8 # Tolerance settings used for odeint. Default values are 1.49012e-8
atol = 1.49012e-8

if choice == 1: # Input variables for the 2-body test case
    G = 6.6748015e-11
    pos = np.array([[-0.5*1.496e11,0.5*1.496e11],[0,0]])
    vel = np.array([[0,0],[-15*1e3,15*1e3]])
    m = np.array([1.989e30, 1.989e30])
    t_tot = 5 * 31557600 # total time to run the simulation for
    dt = 60*60 # time step between outputs calculated by odeint
    animation_tail_length = 200 # number of datapoints to include in the "tail" which trails behind the main body in the animation.
    animation_interval = 20 # delay in milliseconds between animation frames

elif choice == 2: # Input variables for two light bodies and a heavy body - Earth, Mars, Sun
    G = 6.6748015e-11
    pos = np.array([[-1.496e11,-1.524*1.496e11,0],[0,0,0]])
    vel = np.array([[0,0,0],[29.78*1e3,24.07*1e3,0]])
    m = np.array([5.9722e24, 6.4171e23, 1.989e30])
    t_tot = 10 * 31557600 # total time to run the simulation for
    dt = 12*60*60
    animation_tail_length = 80 # number of datapoints to include in the "tail" which trails behind the main body in the animation.
    animation_interval = 40 # delay in milliseconds between animation frames

elif choice == 3: # Input variables for a 3 dimensional four body system
    G = 1
    m = np.array([1000,1,1,1]) # The central body should have a significantly greater mass than the other bodies in the system
    pos = np.array([[0,10,0,0],[0,0,20,0],[0,0,0,30]])
    vel = np.array([[0,0,0,np.sqrt(G*m[0]/pos[2,3])],[0,np.sqrt(G*m[0]/pos[0,1]),0,0],[0,0,np.sqrt(G*m[0]/pos[1,2]),0]]) # For this system, the velocity values are calculated based on the positions so as to give stable orbits
    t_tot = 50 # total time to run the simulation for
    dt = 0.01 # time step between outputs calculated by odeint
    animation_tail_length = 100 # number of datapoints to include in the "tail" which trails behind the main body.
    animation_interval = 2 # delay in milliseconds between animation frames

elif choice in (4,5): # Input variables for Burrau's problem
    G = 1 # Value of the gravitational constant
    pos = np.array([[1,-2,1],[3,-1,-1]]) # initial positions of the three bodies
    vel = np.array([[0,0,0],[0,0,0]]) # initial velocities of the three bodies
    m = np.array([3,4,5]) # masses of the three bodies
    if choice == 5: # Adding an empty set of indices representing the z component
        pos = np.vstack((pos, np.zeros(3)))
        vel = np.vstack((vel, np.zeros(3)))
    t_tot = 70 # total time to run the simulation for
    dt = 0.001 # time step between outputs calculated by odeint
    plot_llim = 40
    plot_ulim = 50
    rtol = 1e-26 # Tolerance settings used for odeint. Default values are 1.49012e-8
    atol = 8.114e-11 # Best value from trial and error
    # atol = 9.549925860215698e-11 # Best value from optimising-burraus-problem.py
    animation_tail_length = 280 # number of datapoints to include in the "tail" which trails behind the main body in the animation.
    animation_interval = 1e-5 # delay in milliseconds between animation frames

else:
    print("Invalid choice, exiting program.")
    exit()

num_dimensions = pos.shape[0]
num_bodies = pos.shape[1]
t_vals = np.arange(0,t_tot+dt,dt) # Array time values to calcuclate the state of the system at.
# For Burrau's problem, we're interested in only plotting specific intervals so that they can be comparted to figures in the paper. These lines enable this.
if choice in (4,5):
    in_range = (t_vals >= plot_llim) & (t_vals <= plot_ulim) # Only points for times between plot_ulim and plot_llim will be plotted.
else:
    in_range = np.arange(t_vals.size) # All points will be plotted

start_time = time.time()
solutions = odeint(solve_func, np.concatenate((pos,vel),axis=None), t_vals, rtol=rtol, atol=atol)
end_time = time.time()
print("Solving the differential equations took",end_time-start_time,"seconds")

fig = plt.figure(1)
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
    if num_dimensions == 2:
        aniax = anifig.add_subplot()
        bodies = aniax.scatter([], [], c=[], marker='.', cmap='tab20')
        anim = FuncAnimation(anifig, animate2d, frames=t_vals.size//10, interval=animation_interval)
    elif num_dimensions == 3:
        aniax = anifig.add_subplot(projection='3d')
        bodies = aniax.scatter([], [], [], c=[], marker='.', cmap='tab20')
        anim = FuncAnimation(anifig, animate3d, frames=t_vals.size//10, interval=animation_interval)
        aniax.set_zlim(1.05*np.min(solutions[:,2*num_bodies:3*num_bodies]), 1.05*np.max(solutions[:,2*num_bodies:3*num_bodies]))
        aniax.set_zlabel("$z$ coordinate")
    else:
        print("Improper number of dimensions, must be 2D or 3D for data to be plotted.")
        exit()
    aniax.set_xlim(1.05*np.min(solutions[:,0:num_bodies]), 1.05*np.max(solutions[:,0:num_bodies]))
    aniax.set_ylim(1.05*np.min(solutions[:,num_bodies:2*num_bodies]), 1.05*np.max(solutions[:,num_bodies:2*num_bodies]))
    aniax.set_xlabel("$x$ coordinate")
    aniax.set_ylabel("$y$ coordinate")
else:
    print("You have chosen not to generate an animation.")

plt.show()
