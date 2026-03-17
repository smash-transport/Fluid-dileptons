import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

datafiles = sys.argv[1:]
lbl = ["pQGP", "lat-QGP FF", "lat-QGP no FF", r"vacuum $\rho$", 'more']
for l,f in zip(lbl,datafiles):
    data = np.loadtxt(f)
    plt.plot(data[:,0],data[:,1],label=l)
plt.yscale('log')
#plt.xlim(0,1.5)
#plt.ylim(5*10**(-9),10**(-4))
#plt.legend(loc=0)
#plt.xlabel("$M$ [GeV]")
#plt.ylabel(r"dR/d$M^2\ [\mathrm{fm}^{-4}\mathrm{GeV}^{-2}]$")
#plt.savefig("latQGP_rate_fig2-1304.2309.pdf")
plt.show()