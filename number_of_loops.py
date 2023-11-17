import numpy as np
import glob
import sys
import matplotlib.pyplot as plt

unloading_rate = ["0.01", "0.001", "0.0001", "0.00001"][::-1]
loading_rate = ["0.01", "0.001", "0.0001", "0.00001"][::-1]

aver_loops = np.zeros((4, 4))
for i, lr in enumerate(loading_rate):
    for j, ur in enumerate(unloading_rate):
        f = glob.glob(sys.argv[1] + f"/{lr}/{ur}/o*500*.o*")
        print(f[0])
        a = np.loadtxt(f[0], skiprows=136, max_rows=5000)
        aver_loops[i, j] = a.mean(axis=0)[-1]


fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(aver_loops - 9999, cmap="viridis")
fig.colorbar(cax)

xaxis = np.arange(len(unloading_rate))
ax.set_xticks(xaxis)
ax.set_yticks(xaxis)
ax.set_xticklabels(loading_rate)
ax.set_yticklabels(unloading_rate)
ax.set_xlabel("loading rate")
ax.set_ylabel("unloading rate")

plt.savefig(sys.argv[1] + "/loop_size.pdf")
