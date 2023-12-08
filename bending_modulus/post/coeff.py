#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

nargv = len(sys.argv)
npoint = 3
if nargv > 2:
    fname = sys.argv[1]
    T = float(sys.argv[2])
    if nargv > 3:
        npoint = int(sys.argv[3])
else:
    print("Usage: coeff.py /path/to/result temperature[K] [# of fitting data]")
    sys.exit()

kT = 1.380649*1e-23*T
data_ = np.genfromtxt(fname)

# sort data
x_ = data_[:,0]
sort_idx = np.argsort(x_)
data = data_[sort_idx]

# calc bending modulus
x = data[:,0]
y = x*x*data[:,1]
c = y[:npoint].mean()
ret = kT/c
print(f"Bending modulus: {ret} [J]")

# plot data
xmin, xmax = 0.02, 0.2
ymin, ymax = 0.0, 0.25
plt.rcParams["font.size"] = 20
plt.rcParams["font.family"] = "Arial"
plt.figure(figsize=(10,7))
plt.scatter(x, y, color="r", s=50)
plt.hlines(c, xmin, xmax, color="gray", linestyle="dashed")
plt.xscale('log')
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.grid(which='major', color='black', linestyle='-')
plt.grid(which='minor', color='gray', linestyle='-')
plt.xlabel(r"$q$ / nm$^{-1}$")
plt.ylabel(r"$q^2<|n^{||}_{q}|^2>$")
plt.tight_layout()
plt.show()
plt.savefig("plot.jpg")
