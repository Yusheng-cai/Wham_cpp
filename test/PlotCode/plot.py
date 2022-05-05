import numpy as np
import matplotlib.pyplot as plt
from scipy.special import logsumexp

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

def DoubleWellPotential(x_range, y_range, numX=50, numY=50, axis=0):
    """
    Function that evaluates the potential at position provided by the user (pos)

    Args:
        xrange(array-like)  : array like object containing [xmin, xmax]
        yrange(array-like)  : array like object containing [ymin, ymax]
        numX(int)           : Number of x in the array (int)
        numY(int)           : Number of y in the array (int)
    """
    x = np.linspace(x_range[0], x_range[1], num=numX)
    y = np.linspace(y_range[0], y_range[1], num=numY)

    xx,yy = np.meshgrid(x,y)
    pos = np.concatenate((xx[:,:,np.newaxis], yy[:,:,np.newaxis]), axis=-1).reshape(-1,2)
    E= 50 * ((pos[:,0] ** 2- 1) **2 + pos[:,1] ** 2).reshape(numX, numY)
    energy1d   = -logsumexp(-E, axis=axis)	

    return x,y, energy1d


if __name__ == "__main__":
    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111)

    fe1 = ReadTabulatedData("../testAdaptive/pjiref.out")
    fe2 = ReadTabulatedData("../testLBFGS/pjiref.out")
    ax.plot(fe1[:,0], fe1[:,3]-fe1[:,3].min(), c='b',label='Adaptive')
    ax.errorbar(fe1[:,0], fe1[:,3]-fe1[:,3].min(),fe1[:,4], c='b', fmt="o", linewidth=3, capsize=10)
    ax.scatter(fe2[:,0], fe2[:,3]-fe2[:,3].min(), c='r', label='LBFGS', s=100)

    x_range = [-1.5,1.5]
    y_range = [-1.5,1.5]
    x,y, DWell = DoubleWellPotential(x_range,y_range)
    DWell -= DWell.min()
    beta = 1000/(8.314 * 300)
    ax.plot(x, beta*DWell, linestyle='dashed', c='black', label='GroundTruth')

    ax.set_xlabel(r"$\tilde{N}$")
    ax.set_ylabel(r"$\beta F$")
    ax.legend(fontsize=15)

    plt.savefig("../Images/validate.png")
