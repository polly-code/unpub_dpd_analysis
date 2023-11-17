import numpy as np
import sys


def msd_pairwise_distance(arr):
    """Mean squared displacement based on pairwise distance, radial MSD"""
    N = len(arr)
    msd = np.zeros(N)
    for i in range(N):
        msd[i] = np.mean(np.square(arr[i + 1 :] - arr[i]))
    return msd


def mean_squared_displacement_1d(arr):
    """Mean squared displacement based on pairwise distance, radial MSD"""
    N = len(arr)
    msd = np.zeros(N)
    for i in range(N):
        msd[i] = np.mean(np.square(arr[i + 1 :] - arr[i]))
    return msd


# parse arguments
d = np.loadtxt(sys.argv[1] + "/dsts.dat", skiprows=100000).T
d = d[80:]
# calculate msd
msds = []
for i in range(len(d)):
    msds.append(msd_pairwise_distance(d[i]))
msds = np.array(msds)
msd = np.mean(msds, axis=0)
# save to file
np.savetxt(sys.argv[1] + "/pairwise_msd.dat", msd)
