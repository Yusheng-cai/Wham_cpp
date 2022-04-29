import numpy as np
import matplotlib.pyplot as plt

from analysis_code.lib.utils import read_dat_gen

if __name__ == "__main__":
    folders = ["testAdaptive", "testLBFGS"]

    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111)

    hist1 = read_dat_gen("testAdaptive/href.out")
    
    N = np.linspace(0,35,hist1.shape[0])
    for i in range(hist1.shape[1]):
        ax.plot(N, hist1[:,i]/hist1[:,i].sum())

    ax.set_xlabel(r"$\tilde{N}$")
    ax.set_ylabel(r"$p(\tilde{N})$")
    ax.legend(fontsize=15)
    plt.savefig("histogram.png")
