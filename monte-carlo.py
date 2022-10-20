import numpy as np
import matplotlib.pyplot as plt

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
    print()

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
