import numpy as np
import matplotlib.pyplot as plt

from analysis_code.lib.utils import read_dat_gen

if __name__ == "__main__":
    folders = ["testAdaptive", "testLBFGS"]

    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111)

    fe1 = read_dat_gen("testAdaptive/pjiref.out")
    fe2 = read_dat_gen("testLBFGS/pjiref.out")
    ax.plot(fe1[:,0], fe1[:,3], c='b',label='Adaptive')
    ax.errorbar(fe1[:,0], fe1[:,3],fe1[:,4], c='b', fmt="o", linewidth=5, capsize=10)
    ax.scatter(fe2[:,0], fe2[:,3], c='r', label='LBFGS', s=60)

    ax.set_xlabel(r"$\tilde{N}$")
    ax.set_ylabel(r"$\beta F$")
    ax.legend(fontsize=15)
    plt.savefig("../Images/validate.png")
