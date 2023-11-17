import scipy
import mdtraj as md
import sys
if 'dcd' in sys.argv[1]:
    tr=md.load_dcd(sys.argv[1], top=sys.argv[2])
else:
    tr=md.load(sys.argv[1], top=sys.argv[2])
tr.center_coordinates()
for i in range(1,len(tr)):
    f1, _ = scipy.spatial.transform.Rotation.align_vectors(tr.xyz[i-1],tr.xyz[i])
    tr.xyz[i] = f1.apply(tr.xyz[i])
tr.save_lammpstrj(sys.argv[1].replace('.lammpstrj','_c.lammpstrj'))
