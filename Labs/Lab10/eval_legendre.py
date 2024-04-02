import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre

def eval_legendre(n,x):
    p = np.zeros(n+1)
    if n == 0:
        p[0] = 1
    elif n == 1:
        p[1] = x
    else:
        p[0] = 1
        p[1] = x
        for i in range(2,n+1):
            p[i] = 1/(i) * ((2*i - 1)*x*p[i-1] - (i -1)*p[i-2])
    return p
    
    
#Example to verify correctness
n = 4
x = 0.75

#My Function
legendre_new = eval_legendre(n,x)

for i in range(5):
    #Scipy Legendre
    legendre_scipy = legendre(i)(x)
    print(f"\ni = {i}\n")
    print("Legendre Polynomial values (My Function):", legendre_new[i])
    print("Legendre Polynomial values (Scipy Function):", legendre_scipy)
    print("\n")