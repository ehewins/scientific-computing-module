import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.special import jn

def task1():
    # For the Bessel functions of this question, we're always working with n=0
    Jn = 1 # Setting up initial value of J0, at z=0
    Js = 0 # Setting up initial value of the slope of J0, at z=0
    dz = 0.05 # Step size for advancing z
    z_start, z_end = 0, 20
    zvals = np.arange(z_start, z_end+dz, dz)
    Jnvals = np.array([Jn])
    for z in zvals[1::]: # Skip the first value of z so we don't divide by zero
        Jnh = Jn + dz/2 * Js
        Jsh = Js + dz/2 * (-Js/z - Jn)
        Jn = Jn + dz*Jsh
        Js = Js + dz*(-Jsh/(z+dz/2) - Jnh)
        Jnvals = np.append(Jnvals, Jn)
    plt.plot(zvals, Jnvals, label='Runga-Kutta method')
    # plt.plot(zvals, integrage.odeint(), label='Scipy odeint function')
    plt.plot(zvals, jn(0,zvals), label='Generated Bessel function')
    plt.legend()
    plt.title("First Bessel function $J_0(z)$")
    plt.xlabel("$z$")
    plt.ylabel("$J_n(z)$")
    plt.show()

def task2():
    g = -9.81 # Accelereation due to gravity in ms^-2.
    initial_speed = 200 # The initial speed of the cannon ball in ms^-1.
    k = 10**-4 # The drag constant in m^-1.
    target_distance = 3000 # Distance to the target in m.
    dt = 0.001 # Size of the timestep for the simulation in s.
    theta = 0 # Starting angle in radians
    dtheta = 0.001 # How much to increase theta by on each iteration.

    target_hit = False
    while target_hit == False and theta <= np.pi/4:
        r_x = 0
        r_y = 0
        v_x = initial_speed * np.cos(theta)
        v_y = initial_speed * np.sin(theta)
        while r_y >= 0:
            r_x += v_x*dt
            r_y += v_y*dt
            v_x += -k*np.sqrt(v_x**2+v_y**2)*v_x*dt
            v_y += (g-k*np.sqrt(v_x**2+v_y**2)*v_y)*dt
        if r_x >= target_distance:
            target_hit = True
        else:
            theta += dtheta

    if target_hit == False:
        print("Could not hit the target, the range was too great.")
    else:
        print("The angle required to hit the target was {:.4f} radians.".format(theta))

def task3():
    # Using units such that ℏ²/2m = ℏ²/2ma² = 1, so we aren't carrying around really small numbers
    V_0 = 10 # The value of the constant V_0, governing the potential inside the well.
    V = lambda x: -V_0 if x <= 1 else 0 # lambda function defining the potential V(x).
    delta_x = 0.001 # The step size between values of position, x, in the calculation
    x_vals = np.arange(0, 7, delta_x) # x values in range 0 <= x < 7, so that any leakage beyond x = a = 1 can be clearly viewed on the graph.
    E = 0 # initial guess for the energy, in terms of ℏ²/2ma².
    delta_E = 0.1 # Starting decrement for each new guess of the energy. The precision of this decrement will be increased as we go.
    plt.figure(1) # Setting up the first figure

    # Finding the ground state energy eigenvalue in 8 refinement stages
    stages = range(1,8+1)
    for stage in stages:
        phi = 0
        phi_vals = np.array([phi]) # We have to set up the φ values array for the initialisation of the following loop.
        while phi_vals[-1] <= 0: # As E is lowered, right as we pass the point of convergence the tail end of the wavefunction goes from below the x axis to above it. We want to stop the refinement stage once we reach this point.
            plt.cla() # Clear the axes from the previous plot
            phi = 0 # Value of φ at x=0, determined by the boundary conditions
            phiprime = 0.01 # Guess for the value of dφ/dx at x=0. It doesn't seem to matter much what you pick for this.
            phi_vals = np.array([phi]) # Set up a fresh φ array, with our initial value already in it.
            for x in x_vals[1::]: # Skipping 1st value of x; we already have the initial value of φ. This loop calculates each next value of φ for the corresponding x.
                phi_h = phi + delta_x/2 * phiprime # Half steps for the O(2) Runge-Kutta method
                phiprime_h = phiprime + delta_x/2 * (V(x)-E)*phi
                phi += delta_x * phiprime_h # Then the full steps
                phiprime += delta_x * (V(x+delta_x)-E)*phi_h
                phi_vals = np.append(phi_vals, phi) # Append this iteration's phi to the array
            plt.plot(x_vals, phi_vals) # Plot the latest version graph
            plt.plot([x_vals[0],x_vals[-1]],[0,0],'k-') # Also plot a line for the x axis
            plt.xlabel("$x$")
            plt.ylabel("$\phi(x)$, unnormalised")
            plt.title("Finding ground state energy: stage {:d}".format(stage))
            plt.pause(0.001)
            print(E) # So we can see what value of E we're on in the console
            E -= delta_E # Next iteration, we try a lower value of E (higher |E|).
        E += 2*delta_E # Reset the value of E to before the last time it was decreased, so that for the next refinement stage we're back on the right side of the x axis.
        delta_E = delta_E / 10 # Increasing the presicion of our E step.
        print("Stage {0:d} complete with E = {1:.8f}".format(stage, E)) # Progress update for the console

    # Finding the area under the graph so that we can normalise it:
    area = integrate.simpson(x_vals, phi_vals)
    phi_vals = phi_vals / np.abs(area)
    plt.figure(2) # Setting up the second figure, for the end results.
    plt.plot(x_vals,phi_vals)
    plt.title("Normalised ground state wavefunction, with $E_0 = ${:.8f}\n(where $\hbar^2/2ma^2 = 1$ and $a=1$)".format(E))
    plt.xlabel("Position, $x$") # Axis labels
    plt.ylabel("Wavefunction value, $\phi(x)$")

    plt.show() # So we actually get to see our graphs.

if __name__ == "__main__":
    choice = int(input("Would you like to run task 1, 2, or 3?\nChoice: "))
    if choice == 1:
        task1()
    elif choice == 2:
        task2()
    elif choice == 3:
        task3()
    else:
        print("Input is not a valid choice.")
