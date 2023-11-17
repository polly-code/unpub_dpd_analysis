import matplotlib.pyplot as plt
import sys
import numpy as np

data = np.loadtxt(sys.argv[1])
if len(sys.argv) == 3:
    mapSize = int(sys.argv[2])
else:
    mapSize = 1000
if len(data) < mapSize:
    print("Size of the data less than desired size of the map")
    sys.exit()
else:
    data = data[:mapSize, :mapSize]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.matshow(
    np.log(data, out=np.zeros_like(data), where=(data != 0)),
    cmap="inferno_r",
    interpolation="none",
)  ## vmin vmax determin the range of colorbar
fig.colorbar(im)
ax.set_xlim(0, mapSize)
ax.set_ylim(mapSize, 0)
fig.savefig(sys.argv[1] + ".pdf")
plt.clf()
