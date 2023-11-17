import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import glob
if sys.argv[1] == '-h':
    print('1st argument is path to files, 2nd argument is starting bead. File mask is /*_c.lammpstrj_msd_s.dat')
    sys.exit()
start_bead = int(sys.argv[2])
path = sys.argv[1]#'/tungstenfs/scratch/ggiorget/pavel/code_projects/plotMSD/1k_test/concentrations/10/'

msds = glob.glob(path+'/*_c.lammpstrj_msd_s.dat')
print(msds)
a=[]
for i in msds:
    a.append(float(i.split('/')[-1].split('_')[0]))
a, msds = zip(*sorted(zip(a, msds)))
plt.figure(figsize=(10,8))
for i in range(len(msds)):
    msd = np.loadtxt(msds[i])
    plt.figure(figsize=(12,10))
    msd=msd[1:].T
    cmap = plt.get_cmap('jet',50)
    col=0
    for numb in range(start_bead,start_bead+50):
        plt.plot(msd[0]*10000, np.mean([msd[j] for j in range(numb,9000,1000)], axis=0),c=cmap(col))
        col += 1
    norm = mpl.colors.Normalize(vmin=start_bead,vmax=start_bead+50)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ticks=np.linspace(start_bead,start_bead+50,11), boundaries=np.arange(start_bead+0.5,start_bead+50.5,1))
    #cbar.set_label('')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('MSDs of 50 beads near CTCF cite')
    plt.show()
    plt.savefig(path+str(a[i])+'_ctcf_msd_50_'+str(start_bead)+'.png')
    
