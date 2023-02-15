# Scientific Computing Third Year Module

This repo is a log of the code I wrote for the scientific computing module I took during the autumn term of the third year of my physics degree.
There are two main components of the module, a series of 5 week-long tasks, and a several week long project (which counts for 80% of the module mark). In my case this project is a Gravitational N-body Modelling

I can all but grantee that when I look back over the code in years to come I'll be amazed by how inefficiently I did everything, but I think it's valuable to have around as a record of my development as a programmer.

## The 5 Week-long Tasks:

### Task 1: Lunar Lander Game

The first task we've been given is a bit of a coding refresh exercise, focused on using matplotlib's widgets module to make a GUI as well as animating a figure.
The idea is simple: you've got a lunar lander falling towards the surface of the moon, which you want to slow down to a safe landing speed. To make it more of a challenge, there should be a limited fuel supply which once depleted disables the thrusters.

### Task 2: Monte Carlo / Stochastic Methods

The second task we've been given covers a series of Stochastic methods, using random numbers to solve physics problems. There are 4 exercises:

1. Generate 10000 uniformly distributed random numbers between 0 and 1 and display them on a histogram. Then transform the random numbers into a new set with a linearly increasing probability density function between 0 and 1 which should then also be shown on a histogram.

2. Use Monte Carlo integration methods to evaluate the integral of exp(-|x|)cos(x)dx between the limits x=-π/2 and x=π/2. Then compare the result to one calculated with an inbuilt integration function.

3. Use a Monte Carlo method to calculate the volume of a torus with outer radius 10cm and inner radius 5cm. The equation for the surface of a torus with outer radius b and inner radius a is: (√(x²+y²) - (b+a)/2)² + z² = ((b-a)/2)²

4. Write a program to simulate a random walker in two-dimensions. At each time step the walker can move in one of the four possible directions (up, down, left, right) with equal probability. By considering the trajectories of many walkers, determine how the root mean square distance from the origin increases with time.

### Task 3: Solving Differential Equations / Ordinary Differential Equations

The third task we've been set is to solve 3 ordinary differential equation problems using numerical methods. The third of these three sub-tasks is assessed content, so I'll be putting extra work into that one. The 3 exercises are as follows:

1. Solve Bessel's equation using the 2nd order Runga-Kutta method for n=0 and the boundary conditions J₀(0) = 1 and dJ₀/dz |_{z=0} = 0.
   Bessel's equation: z² * (d²Jₙ)/(dz²) + z * (dJₙ)/(dz) + (z²-n²)*Jₙ = 0.
   Then find the solution using an inbuilt ODE solver and compare your solutions with one from the inbuilt Bessel function generator.

2. Use the shooting method to find the angle of elevation a cannon must be fired at for a cannon ball with initial speed 200ms⁻¹ to hit a target 3km away, assuming flat ground. The equation of motion is:
   m (d²r)/(dt²) = mg - mk|v|v
   where r is the position vector, v is the velocity vector and g is the acceleration due to gravity. The motion may be assumed to be two-dimensional and |v| = √(v_x² = v_y²) is the magnitude of the velocity and the drag constant k is 10⁻⁴ m⁻¹.

3. The potential energy function for a particle confined by a surface can be written in the form:
           /              ∞,  x <= 0
   V(x) = {  -(ℏ²V₀)/(2ma²),  0 < x <= a
           \              0,  x > a
   Where V₀ is a constant. The wavefunction satisfies the time-independent Schrödinger equation (-ℏ²/2m)(d²φ/dx²) + V(x)φ = Eφ.
   If V₀ = 10, find the ground state energy (in terms of ℏ²/2ma²) and normalised ground state wavefunction using the shooting method. As the potential is infinite for x<=0, the wavefunction will be zero for all negative x.

The mark received for Task 3 was 5 out of 6 marks. No further feedback has been given unfortunately.

### Task 4: Signal Processing

The fourth task is all about signal processing using Fourier methods. There's also a focus on GUIs again, which we'll be using to change various parameters such as the signal frequency and sampling rate so their effect on the power spectrum can be observed in real time. The GUI should generate a signal and plot that signal against time. First introduce a slider that controls the frequency of a sinusoidal signal, and a close button which closes the GUI. There should also be labels on the figure and slider to indicate the units. Calculate the Fourier transform of the signal, and use this to plot a power spectrum that updates when slider is used to change the frequency.

The following are suggestions for additional functionality which can be added once the base program is working:

* Adding sliders controlling, for example, the number of sampling points, total time, amplitude and phase of the signal
* Adding additional signal types (square wave, saw tooth)
* Windowing in the time domain
* Adding noise
* Adding multiple signals
* Windowing in the frequency space, used for plotting a "cleaned-up" inverse Fourier transform
* Filtering / convolution
* 2D Fourier transforms
* Data compression

NOTE: This week I had a lot of my time taken up by one of my other modules, so I didn't spend as much time on this one as I normally like to. As a result, I hadn't managed to add all the features I was hoping to by the time the deadline came around, so the code in signal-processing.py still has some fearsome bugs when it comes to the windowing (the bits where you can change the upper and lower limits of the time and frequency) and the inverse Fourier transform (which only behaves itself when you don't window the data at all, i.e. when the output is identical to the original signal). If the lower limit and upper limit are set to the very ends of their ranges, the program works fine however. I was considering removing these bits of code, but I think I'll leave them in and aim to go back and fix it at some later stage. Not sure when I'll get time to, but I hope I can.

The mark received for Task 4 was 6 out of 7 marks. No further feedback was given.

### Task 5: Optimisation

This week's tasks are themed around optimising code, and learning best practices for writing fast. We'll be given a badly written program (called neighbour.py) which we're tasked with improving by only making changes between the lines of code beginning with "start_time" and "end_time". The goal is to optimise the code to the point where it runs in less than 12 seconds for N=20000 and seed=1234 on the computers in the lab. (If I can get it to execute in less than 12 seconds on this old laptop, then it's probably safe to assume it'll run faster than that on the desktops in the lab). The file `optimisation.py` will start out as a copy of neighbour.py, and you should be able to track the changes I make to it through the git log.

The mark received for Task 5 was 7 out of 7 marks. Full marks!

## Project: Gravitational N-body Modelling

The project I've chosen to work on is the N-body simulation, for which I shall solve by direct summation the gravitational acceleration between each body/particle i and all other bodies/particles j, in order to determine the evolution of the N-body system through time.

(d²rᵢ)/(dt²) = ∑ᴺⱼ₌₁ (Gmⱼrᵢⱼ)/(|rᵢⱼ|³)

Where rᵢ is the position vector of the object i, rᵢⱼ is the vector from i's position to j's position, i.e. rᵢⱼ = rⱼ-rᵢ, mⱼ is the mass of object j, G is the gravitational constant, t is time and N is the total number of bodies/objects/particles in the simulation.
This is a second order differential equation, which will need to be split into a pair of coupled first order differential equations before they are computed by an appropriate ODE solving python module.

The directory N_body_simulation contains all the files I submitted for grading. The README and pdf of figures were requirements for the submission.

### Project Progression Guidelines:

1. Begin by simulating 2 bodies in 2 dimensions and check that the results lead to simple, closed, elliptical orbits.
2. Check that energy is conserved by the simulation. With the default tolerance settings of scipy's odeint solver, you should find that |ΔE| = (|E(t)-E(0)|)/(E(0)) < 10⁻⁵.
3. Extend this 2D code to simulate three bodies, and compare the results with some of the following analytically solvable cases:
    * Two light bodies in orbit around a much more massive one. Should reproduce simple elliptical orbits.
    * Two heavy bodies in close orbit with a light body orbiting them both. Can investigate the behaviour of the light body as a function of its position relative to the heavy bodies, as well as identify the Lagrangian points.
    * Solve Burrau's problem, see the paper "Complete Solution of a General Problem of Three Bodies" by Szebehely and Peters, DOI: 10.1086/110355 .
4. Extend the code to three dimensions and many bodies. Try out some of the following more advanced simulations:
    * Study the stability of known exoplanetary systems, looking up the relevant planetary data to do so.
    * Study the stability of the solar system under close encounters with external objects (such as rogue stars or planets passing through the solar system), or the effect of an additional large planet.
    * Study the evolution of a star cluster.
