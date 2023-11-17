import matplotlib.pyplot as plt
import sys
import numpy as np

data = np.loadtxt(sys.argv[1])
mapSize = int(sys.argv[2])
data = data[:mapSize, :mapSize]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.matshow(
    np.log(data, out=np.zeros_like(data), where=(data != 0)),
    cmap="inferno",
    interpolation="none",
)  ## vmin vmax determin the range of colorbar
fig.colorbar(im)
# ax.plot([450,450],[450,1000], '--', c='b')
# ax.plot([550,550],[550,1000], '--', c='b')

ax.set_xlim(0, mapSize)
ax.set_ylim(mapSize, 0)
# ax.set_xbound(lower=1000, upper=0)
# ax.set_ybound(lower=1000, upper=0)
fig.savefig(sys.argv[1] + ".png", dpi=300)
plt.clf()
