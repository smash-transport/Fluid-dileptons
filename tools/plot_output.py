import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

datafiles = sys.argv[1:]
labels = {9001: r"$\rho$", 9002: r"$\omega$", 9003: r"$\phi$", 9004: r"multi-$\pi$", 9005: r"QGP"}
sum_data = {}

for f in datafiles:
    data = np.loadtxt(f)
    unique_vals = labels.keys()
    for val in unique_vals:
        mask = data[:,2] == val
        filtered = data[mask]
        plt.plot(filtered[:,0], filtered[:,1], label=labels[val])
        # Accumulate sum for each x-value
        for x, y in zip(filtered[:,0], filtered[:,1]):
            x_key = float(x)
            if x_key not in sum_data:
                sum_data[x_key] = 0
            sum_data[x_key] += y

if sum_data:
    x_vals = sorted(sum_data.keys())
    y_vals = [sum_data[x] for x in x_vals]
    plt.plot(x_vals, y_vals, color='black', linewidth=2, label='sum')

plt.yscale('log')
#plt.xlim(0,1.5)
plt.ylim(5*10**(-17),10**(-4))
plt.legend(loc=0)
plt.xlabel("$M$ [GeV]")
plt.ylabel(r"dR/d$M\ [\mathrm{fm}^{-4}\mathrm{GeV^{-1}}]$")
#plt.savefig("latQGP_rate_fig2-1304.2309.pdf")
plt.show()