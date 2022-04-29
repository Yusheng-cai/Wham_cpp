import numpy as np
import matplotlib.pyplot as plt

from analysis_code.lib.utils import read_dat_gen

if __name__ == "__main__":
    folders = ["testAdaptive", "testLBFGS"]

    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111)

    kl1 = read_dat_gen("testAdaptive/klref.out")
    x = np.arange(1,len(kl1)+1)

    ax.bar(x, kl1.flatten())
    max_ =kl1.max() * 4
    ax.set_ylim(0,max_)

    ax.set_xlabel(r"Simulation")
    ax.set_ylabel(r"KL Divergence")
    plt.savefig("../Images/KL.png")
