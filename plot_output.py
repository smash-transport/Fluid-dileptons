import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

datafiles = sys.argv[1:]
lbl = ["pQGP", "lat-QGP FF", "lat-QGP no FF", r"vacuum $\rho$", 'more']
sum_data = {}  # Dictionary to store sum of y-values for each x-value

for l,f in zip(lbl,datafiles):
    data = np.loadtxt(f)
    unique_vals = np.unique(data[:,2])
    for val in unique_vals:
        mask = data[:,2] == val
        filtered_data = data[mask]
        # Plot with label indicating the value from column 2
        plt.plot(filtered_data[:,0], filtered_data[:,1], label=f'{l} (col3={val})')
        # Accumulate sum for each x-value
        for x, y in zip(filtered_data[:,0], filtered_data[:,1]):
            x_key = float(x)
            if x_key not in sum_data:
                sum_data[x_key] = 0
            sum_data[x_key] += y

# Plot the sum of all curves in black
if sum_data:
    x_vals = sorted(sum_data.keys())
    y_vals = [sum_data[x] for x in x_vals]
    plt.plot(x_vals, y_vals, color='black', linewidth=2, label='Sum')

plt.yscale('log')
#plt.xlim(0,1.5)
plt.ylim(5*10**(-17),10**(-4))
#plt.legend(loc=0)
#plt.xlabel("$M$ [GeV]")
#plt.ylabel(r"dR/d$M^2\ [\mathrm{fm}^{-4}\mathrm{GeV}^{-2}]$")
#plt.savefig("latQGP_rate_fig2-1304.2309.pdf")
plt.show()