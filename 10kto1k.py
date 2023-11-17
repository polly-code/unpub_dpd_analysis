import sys
import numpy as np
a=np.loadtxt(sys.argv[1])
if len(sys.argv) == 3:
    aver = True
else:
    aver = False
init = True
for i in range(1,9):
    if init:
        b = a[i*1000:(i+1)*1000, i*1000:(i+1)*1000]
        init = False
    else:
        b += a[i*1000:(i+1)*1000, i*1000:(i+1)*1000]
if aver:
    np.savetxt(sys.argv[1]+'_s', b/8)
else:
    np.savetxt(sys.argv[1]+'_s', b)