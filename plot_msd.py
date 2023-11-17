import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

path = sys.argv[1]
msds = glob.glob(path + "/*_c.lammpstrj0_msd_a.dat")
a = []
for i in msds:
    a.append(float(i.split("/")[-1].split("_")[0]))
print(a)
sys.exit()
a, msds = zip(*sorted(zip(a, msds)))
for i in range(len(a)):
    print(a[i], msds[i])
plt.figure(figsize=(10, 8))
for i in range(len(msds)):
    msd = np.loadtxt(msds[i])
    plt.plot(msd[:, 0] * 10000, msd[:, 1], label="le " + str(a[i]))
    x = np.log10(msd[0:5, 0] * 10000)
    y = np.log10(msd[0:5, 1])
    slope, lgd = np.polyfit(x, y, 1)
    plt.gcf().text(
        x=0.15,
        y=0.55 - i / 40,
        s=r"$\alpha_{"
        + str(a[i])
        + "} = $"
        + "{:.3f}".format(slope)
        + r" D$_{"
        + str(a[i])
        + "}$ = "
        + "{:.4f}".format(10 ** lgd)
        + r" - $10^4$-$5\cdot10^4$ a.u.",
    )
    x = np.log10(msd[50:100, 0] * 10000)
    y = np.log10(msd[50:100, 1])
    slope, lgd = np.polyfit(x, y, 1)
    plt.gcf().text(
        x=0.15,
        y=0.85 - i / 40,
        s=r"$\alpha_{"
        + str(a[i])
        + "} = $"
        + "{:.3f}".format(slope)
        + r" D$_{"
        + str(a[i])
        + "}$ = "
        + "{:.4f}".format(10 ** lgd)
        + " - $5\cdot10^5$-$10^6$ a.u.",
    )

# path='concentrations/'
# msds = glob.glob(path+'/*_c.lammpstrj0_msd_a.dat')
# msds = np.sort(msds)
# print(msds)

# for i in [1]:#range(len(msds)):
#    msd = np.loadtxt(msds[i])
#    print(msd[:,0])
#    plt.plot(msd[:,0]*10000, msd[:,1], '--', label=str(i))
#    x = np.log10(msd[0:5,0]*10000)
#    y = np.log10(msd[0:5,1])
#    slope, lgd = np.polyfit(x, y, 1)
#    plt.gcf().text(x= 0.15, y= 0.7, s = r'$\alpha_{NCh} = $' + '{:.3f}'.format(slope)+r' D$_{NCh}$ = ' + '{:.4f}'.format(10**lgd) + ' - 10-50 a.u.')
#    x = np.log10(msd[50:100,0]*10000)
#    y = np.log10(msd[50:100,1])
#    slope, lgd = np.polyfit(x, y, 1)
#    plt.gcf().text(x= 0.15, y= 0.75, s = r'$\alpha_{LE} = $' + '{:.3f}'.format(slope)+r' D$_{LE}$ = ' + '{:.4f}'.format(10**lgd) + ' - 10-50 a.u.')

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\Delta$t, a.u.")
plt.ylabel(r"MSD, $a.u.^2$")
plt.legend(loc="lower right")
plt.savefig(sys.argv[1] + "msd_aver.png", dpi=300)
