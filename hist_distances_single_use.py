import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys
import pandas as pd
import scipy.stats as ss
import itertools

f = open(sys.argv[1] + "dsts.dat", "r")
for i in f:
    header = i.split()
    print("header is", header, len(header))
    f.close()
    break
d_on = np.loadtxt(sys.argv[1] + "dsts.dat", skiprows=100000)  # , max_rows=1000)
d_on = d_on.T
print("length of dsts ", len(d_on))

f = open((sys.argv[1] + "dsts.dat").replace("/on/", "/off/"), "r")
for i in f:
    header = i.split()
    print("header is", header, len(header))
    f.close()
    break
d_off = np.loadtxt(
    (sys.argv[1] + "dsts.dat").replace("/on/", "/off/"), skiprows=100000
)  # , max_rows=1000)
d_off = d_off.T
print("length of dsts ", len(d_off))


ctcf_on = []
for i in range(80, 88):
    ctcf_on.append(d_on[i])
    print(header[i])
plt.figure(figsize=(8, 6))
hson = []
for i in range(8):
    h, b = np.histogram(ctcf_on[i], bins=50, range=(0, 12))
    hson.append(h)
ctcf_off = []
for i in range(80, 88):
    ctcf_off.append(d_off[i])
    print(header[i])
plt.figure(figsize=(8, 6))
hsoff = []
for i in range(8):
    h, b = np.histogram(ctcf_off[i], bins=50, range=(0, 12))
    hsoff.append(h)
plt.plot((b[:-1] + b[1:]) / 2, np.mean(hson, axis=0), label="ctcf on")
plt.plot((b[:-1] + b[1:]) / 2, np.mean(hsoff, axis=0), label="ctcf off")
plt.title("Between CTCFs")
plt.xlabel("Distance, a.u.")
plt.ylabel("Count")
plt.legend()
plt.savefig(sys.argv[1] + "hist_ctcf_on_off.png", dpi=300)
plt.clf()

sys.exit()
