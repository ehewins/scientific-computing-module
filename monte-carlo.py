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
    print("With a ransom sample of size {0:d}, the integral evaluated using Monte Carlo methods gives\nI = {1:.4f}, and when calculated numerically, gives I = {2:.4f}.".format(N, I, trueI))

def task3():
    print()

def task4():
    print()

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
