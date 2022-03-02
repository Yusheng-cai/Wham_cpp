import numpy as np
import matplotlib.pyplot as plt

from analysis_code.lib.utils import read_dat_gen

if __name__ == "__main__":
    folders = ["testAdaptive", "testLBFGS"]

    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111)



    fe1 = read_dat_gen("testAdaptive/pjiref.out")
    fe2 = read_dat_gen("testLBFGS/pjiref.out")
    fe3 = read_dat_gen("testSparseSampling/pjiref.out")
    ax.plot(fe1[:,0], fe1[:,-1], label="Adaptive", c='b')
    ax.scatter(fe2[:,0], fe2[:,-1], label='LBFGS', c='r', s=60)
    ax.scatter(fe3[:,0], fe3[:,-1], label='Sparse', c='g', s=60)

    ax.set_xlabel(r"$\tilde{N}$")
    ax.set_ylabel(r"$\beta F$")
    ax.legend(fontsize=15)
    plt.savefig("validate.png")
