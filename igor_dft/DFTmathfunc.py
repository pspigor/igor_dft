import numpy as np

def MATH_TrapezoidalRule(f,a,b) :
    N = np.size(f)
    h = (b-a)/(N-1)
    W = h*np.ones(N)
    W[ 0] = h/2
    W[-1] = h/2
    I = np.sum(W*f)
    return I
