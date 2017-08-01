#!/usr/bin/env python
import numpy as np
import math

def generate():
    f = open('data.dat', 'w')
    X = np.random.uniform(-1.0, 1.0, 100)
    Y = np.random.uniform(-1.0, 1.0, 100)
    value = np.exp(-(X**2 + Y**2))
    value += np.random.normal(0.0, 0.2, 100)
    error = np.random.uniform(0.1, 0.2, 100) * value
    f.write('x\t y\t value\t error\n')
    for i in range(100):
        f.write('%.6f\t %.6f\t %.6f\t %.6f\n' %(X[i], Y[i], value[i], error[i]))
    f.close()
    return

if __name__ == '__main__':
    generate()

