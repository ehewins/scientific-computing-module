# Scientific Computing Third Year Module

This repo is a log of the code I wrote for the scientific computing module I took during the autumn term of the third year of my physics degree.

I can all but grantee that when I look back over the code in years to come I'll be amazed by how inefficiently I did everything, but I think it's valuable to have around as a record of my development as a programmer.

## First Task: Lunar Lander Game

The first task we've been given is a bit of a coding refresh exercise, focused on using matplotlib's widgets module to make a GUI as well as animating a figure.
The idea is simple: you've got a lunar lander falling towards the surface of the moon, which you want to slow down to a safe landing speed. To make it more of a challenge, there should be a limited fuel supply which once depleted disables the thrusters.

## Second Task: Monte Carlo / Stochastic Methods

The second task we've been given covers a series of Stochastic methods, using random numbers to solve physics problems. There are 4 exercises:

1. Generate 10000 uniformly distributed random numbers between 0 and 1 and display them on a histogram. Then transform the random numbers into a new set with a linearly increasing probability density function between 0 and 1 which should then also be shown on a histogram.
2. Use Monte Carlo integration methods to evaluate the integral of exp(-|x|)cos(x)dx between the limits x=-π/2 and x=π/2. Then compare the result to one calculated with an inbuilt integration function.
3. Use a Monte Carlo method to calculate the volume of a torus with outer radius 10cm and inner radius 5cm. The equation for the surface of a torus with outer radius b and inner radius a is: ((x^2+y^2)^0.5 - (b+a)/2)^2 + z^2 = ((b-a)/2)^2
4. Write a program to simulate a random walker in two-dimensions. At each time step the walker can move in one of the four possible directions (up, down, left, right) with equal probability. By considering the trajectories of many walkers, determine how the root mean square distance from the origin increases with time.
