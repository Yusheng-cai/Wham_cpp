import numpy as np
import matplotlib.pyplot as plt

def ReadTabulatedData(file_path):
    """
    Function that reads the .dat from INDUS simulations

    Args:
        file_path(str): the path to the file 

    Return:
        an numpy array that contains the N and Ntilde from the INDUS simulation (N, Ntilde) -> where both are of shape (nobs,)
    """
    f = open(file_path)
    lines = f.readlines()
    lines = [line for line in lines if line[0]!="#"]

    lines = [[float(num) for num in line.rstrip("\n").split()] for line in lines if line[0]!="#"]
    lines = np.array(lines)

    return lines



if __name__ == "__main__":
    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111)

    kl1 = ReadTabulatedData("../testAdaptive/klref.out")
    x = np.arange(1,len(kl1)+1)

    ax.bar(x, np.abs(kl1.flatten()))
    max_ =kl1.max() * 4
    ax.set_ylim(0,max_)

    ax.set_xlabel(r"Simulation")
    ax.set_ylabel(r"KL Divergence")
    plt.savefig("../Images/KL.png")
