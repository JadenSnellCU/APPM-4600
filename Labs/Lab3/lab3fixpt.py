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
    
x = 7 ** (1/5)
# use routines 
f1 = lambda x: x*(1+ ((7- x**5)/x**2))**3



f2 = lambda x: x-((x**5 - 7)/x**2)


f3 = lambda x: x-((x**5 - 7)/5*x**4)


f4 = lambda x: x-((x**5 - 7)/12)

Nmax = 100
tol = 10e-10

#Problems a-c commit overflow
'''
x0 = 1
[xstar,ier] = fixedpt(f1,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f1(xstar):',f1(xstar))
print('Error message reads:',ier)


x0 = 1
[xstar,ier] = fixedpt(f2,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f2(xstar):',f2(xstar))
print('Error message reads:',ier)


x0 = 1
[xstar,ier] = fixedpt(f3,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f3(xstar):',f3(xstar))
print('Error message reads:',ier)

'''

x0 = 1
[xstar,ier] = fixedpt(f4,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f4(xstar):',f4(xstar))
print('Error message reads:',ier)

