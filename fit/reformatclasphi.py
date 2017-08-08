#!/usr/bin/env python
import numpy as np
import os, sys

def Reformat(fs):
    data = np.loadtxt(fs, skiprows = 8)
    f = open(fs.replace(".txt", ".dat"), "w")
    data[:,1] = ((data[:,1] + 6e-3) // 0.025) * 0.025
    for d in data:
        f.write("%.3f\t%.3f\t%.6f\t%.6f\n" %(d[0], d[1], d[2], d[3]))
    f.close()
    return

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit()
    else:
        Reformat(sys.argv[1])
        sys.exit()
