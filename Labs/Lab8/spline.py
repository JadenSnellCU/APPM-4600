import numpy as np
import matplotlib.pyplot as plt
import math
from numpy.linalg import inv

def eval_line(x0,x1,fx0,fx1,alpha):
    slope = (fx1-fx0)/(x1-x0)
    y_a = fx0 + slope * (alpha - x0)
    
    return y_a

def eval_lin_spline(xeval,Neval,a,b,f,Nint):
    xint = np.linspace(a,b,Nint +1)
    yeval = np.zeros(Neval)
    for j in range(Nint):
        atmp = xint[j]
        btmp = xint[j+1]
        
        ind = np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]
        n = len(xloc)
        fa = f(atmp)
        fb = f(btmp)
        yloc = np.zeros(len(xloc))
        
        for kk in range(n):
            yloc[kk] = eval_line(atmp,btmp,fa,fb,xloc[kk])
            
            
        yeval[ind] = yloc
        return yeval

f = lambda x : 1 / (1+(10*x)**2)

a = -1
b = 1

Neval = 100
xeval = np.linspace(a,b,Neval)

Nint = 10

yeval = eval_lin_spline(xeval,Neval,a,b,f,Nint)

fex = f(xeval)

plt.figure()
plt.plot(xeval, fex, 'ro-', label='Original Function')
plt.plot(xeval, yeval, 'bs-', label='Spline')
plt.legend()
plt.show()

err = abs(yeval - fex)
plt.figure()
plt.plot(xeval, err, 'ro-')
plt.title('Error Plot')
plt.show()

