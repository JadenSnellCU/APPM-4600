# import libraries
import numpy as np
    
# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]
    

# use routines 
f1 = lambda x: -np.sin(2*x) + 1.25*x - 0.75
''' 
fixed point is alpha1 = 1.4987....
'''

Nmax = 100
tol = 1e-10

''' test f1 '''
    
x0 = 0
print('\n')
[xstar,ier] = fixedpt(f1,x0,tol,Nmax)
print('Root 1:',round(xstar,10))
print('Error message reads:',ier)
print('\n')

x0 = 3
[xstar,ier] = fixedpt(f1,x0,tol,Nmax)
print('Root 2:',round(xstar,10))
print('Error message reads:',ier)
print('\n')

