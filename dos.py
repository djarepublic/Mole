import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def DOS(file, atom):
    f = open(file, 'r')
    lines = f.readlines()

    line = 'ion ' + str(atom) 
    st_mas = []
    for i in range(len(lines)):
        if line in lines[i]:
            print(i)
            i = i + 2
            break

    for j in range(i, 10**7):
        if 'set' in lines[j]:
            break
        numbers = re.findall(r'-?\d+(?:\.\d+)?', lines[j])
        numbers = [float(num) for num in numbers]
        st_mas.append(numbers)

    mas = np.array(st_mas)
    en  = np.array([row[0] for row in mas])
    s = np.array([row[1] for row in mas])
    py = np.array([row[2] for row in mas])
    pz = np.array([row[3] for row in mas])
    px = np.array([row[4] for row in mas])
    dxy = np.array([row[5] for row in mas])
    dyz = np.array([row[6] for row in mas])
    dz2 = np.array([row[7] for row in mas])
    dxz = np.array([row[8] for row in mas])
    dx2y2 = np.array([row[9] for row in mas])
    fy3x2 = np.array([row[10] for row in mas])
    fxyz = np.array([row[11] for row in mas])
    fyz2 = np.array([row[12] for row in mas])
    fz3 = np.array([row[13] for row in mas])
    fxz2 = np.array([row[14] for row in mas])
    fzx2 = np.array([row[15] for row in mas])
    fx3 = np.array([row[16] for row in mas])
    p = py + pz + px
    d = dxy + dyz + dz2 + dxz + dx2y2
    f = fy3x2 + fxyz + fyz2 + fz3 + fxz2 + fzx2 + fx3
    tot = f + d + p

    orb0 = [s, p, d, f, tot]
    orb = [en]
    for o in orb0:
        oi = savgol_filter(o, 17, 1)
        orb.append(oi)
    
    return(orb)