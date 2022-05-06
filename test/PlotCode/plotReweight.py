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

    reweightAvg = ReadTabulatedData("../testReweight/ReweightAvgRef.out")
    beta = 1000/(8.314 * 300)
    x = beta * np.linspace(-100,100,1000)
    
    ax.plot(x, reweightAvg[:,1], c='b', linewidth=5)
    ax.set_xlabel(r"$\beta \phi$")
    ax.set_ylabel(r"$\langle x \rangle$")
    plt.savefig("../Images/Reweight.png")
