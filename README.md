# Scientific Computing Third Year Module

This repo is a log of the code I wrote for the scientific computing module I took during the autumn term of the third year of my physics degree.

I can all but grantee that when I look back over the code in years to come I'll be amazed by how inefficiently I did everything, but I think it's valuable to have around as a record of my development as a programmer.

## Task 1: Lunar Lander Game

The first task we've been given is a bit of a coding refresh exercise, focused on using matplotlib's widgets module to make a GUI as well as animating a figure.
The idea is simple: you've got a lunar lander falling towards the surface of the moon, which you want to slow down to a safe landing speed. To make it more of a challenge, there should be a limited fuel supply which once depleted disables the thrusters.

## Task 2: Monte Carlo / Stochastic Methods

The second task we've been given covers a series of Stochastic methods, using random numbers to solve physics problems. There are 4 exercises:

1. Generate 10000 uniformly distributed random numbers between 0 and 1 and display them on a histogram. Then transform the random numbers into a new set with a linearly increasing probability density function between 0 and 1 which should then also be shown on a histogram.

2. Use Monte Carlo integration methods to evaluate the integral of exp(-|x|)cos(x)dx between the limits x=-π/2 and x=π/2. Then compare the result to one calculated with an inbuilt integration function.

3. Use a Monte Carlo method to calculate the volume of a torus with outer radius 10cm and inner radius 5cm. The equation for the surface of a torus with outer radius b and inner radius a is: (√(x²+y²) - (b+a)/2)² + z² = ((b-a)/2)²

4. Write a program to simulate a random walker in two-dimensions. At each time step the walker can move in one of the four possible directions (up, down, left, right) with equal probability. By considering the trajectories of many walkers, determine how the root mean square distance from the origin increases with time.

## Task 3: Solving Differential Equations / Ordinary Differential Equations

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
