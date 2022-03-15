import matplotlib.pyplot as plt
from wham.lib.wham_utils import read_dat_gen
import numpy as np

if __name__ == "__main__":
    data = read_dat_gen("a.out")
    ax   = plt.figure().add_subplot(111)
    phi  = np.linspace(-10,10,500)

    ax.plot(data[:,0], data[:,1])
    plt.savefig("test.png")
