import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def task1():
    uniform_rand = np.random.random_sample(10000)
    plt.hist(uniform_rand, bins=20)
    plt.title("Histogram of 10000 uniformly distributed numbers [0.,1.)")
    plt.figure()
    linear_rand = np.sqrt(uniform_rand)
    plt.hist(linear_rand,bins=20)
    plt.title("Histogram of 10000 numbers\nwith a linearly increasing prob. density [0.,1.)")
    plt.show()

def task2():
    f = lambda x: np.exp(-np.abs(x))*np.cos(x)
    lowerlim, upperlim = -np.pi/2, np.pi/2
    N = 10000
    V = upperlim - lowerlim
    uniform_rand = lowerlim + V*np.random.random_sample(N)
    I = V/N * np.sum(f(uniform_rand))
    trueI,_ = integrate.quad(f, lowerlim, upperlim)
    percentdiff = 100*np.abs(I-trueI)/trueI
    print("With a ransom sample of size {0:d}, the integral evaluated using Monte Carlo methods gives\nI = {1:.4f}, and when calculated numerically, gives I = {2:.4f}.\nThese results differ by {2:.3f}%.".format(N, I, trueI, percentdiff))

def task3():
    a, b = 5, 10 # Inner and outer radii of the torus in cm
    # Generate a random collection of 3D coordinates, within an appropriate sized volume
    N = 100000
    xcoords = -b + 2*b*np.random.random_sample(N)
    ycoords = -b + 2*b*np.random.random_sample(N)
    zcoords = -(b-a)/2 + (b-a)*np.random.random_sample(N)
    incount = 0
    for i in range(N):
        x, y, z = xcoords[i], ycoords[i], zcoords[i]
        if (np.sqrt(x**2+y**2) - (b+a)/2)**2 + z**2 < ((b-a)/2)**2:
            incount += 1
    boxV = 2*b * 2*b * (b-a)
    torusV = boxV * incount/N
    print("The volume of a torus with inner radius {0:.1f}cm and outer radius {1:.1f}cm was found to be approximately {2:.2f}cm^3 by Monte Carlo methods.".format(a,b,torusV))
    calcV = np.pi**2/4*(b+a)*(b-a)**2
    print("The true volume is {:.2f}cm^3.".format(calcV))

def task4():
    walks = 1000 # Number of random walks averaged over
    steps = 500 # Number of steps per walk
    xhist = np.zeros(walks)
    yhist = np.zeros(walks)
    rms_rhist = np.array([])
    for i in range(steps):
        for w in range(walks):
            x = xhist[w]
            y = yhist[w]
            match np.random.choice(["up", "down", "left", "right"]):
                case "up":
                    y += 1
                case "down":
                    y -= 1
                case "left":
                    x -= 1
                case "right":
                    x += 1
            xhist[w] = x
            yhist[w] = y
        rms_r = np.sqrt(np.mean(xhist**2 + yhist**2))
        rms_rhist = np.append(rms_rhist, rms_r)
        print("Running: {:.1f}% complete".format(100*i/steps))

    plt.title("RMS distance from origin of {:d} random walks.".format(walks))
    plt.xlabel("Step number")
    plt.ylabel("RMS distance from origin")
    plt.plot(range(steps)+1, rms_rhist)
    plt.show()


if __name__ == "__main__":
    choice = int(input("Would you like to run task 1, 2, 3, or 4?\nChoice: "))
    if choice == 1:
        task1()
    elif choice == 2:
        task2()
    elif choice == 3:
        task3()
    elif choice == 4:
        task4()
    else:
        print("Input is not a valid choice.")
