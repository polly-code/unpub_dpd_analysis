import numpy as np
import os

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

ctcf = ['on', 'off']
extruder_speed = ['40', '200', '500', '1000']
loading_prob = ['0.01', '0.001', '0.0001', '0.00001']
unloading_prob = ['0.01', '0.001', '0.0001', '0.00001']
total_sys = []
for a in ctcf:
    for b in extruder_speed:
        for c in loading_prob:
            for d in unloading_prob:
                total_sys.append([a,b,c,d])
                
for i in range(len(total_sys)):
    path = 'data/' + '/'.join(total_sys[i])
    
    check = any(os.scandir(path))
    total_sys[i].append(check)
    if check:
        print(path, color.GREEN + color.BOLD + 'True' + color.END)
    else:
        print(path, color.RED + color.BOLD + 'False' + color.END)
