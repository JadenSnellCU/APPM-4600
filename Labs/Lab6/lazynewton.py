import numpy as np
import math
import time
from numpy.linalg import inv
from numpy.linalg import norm
import matplotlib.pyplot as plt;



def evalF(x):

    F = np.zeros(2)

    F[0] = 4*x[0]**2+x[1]**2 - 4
    F[1] = x[0]+x[1] - np.sin(x[0] - x[1])
    return F

def evalJ(x):
    J = np.array([[8*x[0], 2*x[1]],[1-np.cos(x[0]-x[1]), 1 + np.cos(x[0]-x[1])]])
    return J

def LazyNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    xlist = np.zeros((Nmax+1,len(x0)))
    xlist[0] = x0
    J = evalJ(x0)
    
    for its in range(Nmax):
       if its % 2 == 0:
           J = evalJ(x0) 
       F = evalF(x0)
       x1 = x0 - np.linalg.solve(J,F)
       xlist[its+1]=x1

       if (norm(x1-x0) < tol*norm(x0)):
           xstar = x1
           ier =0
           return[xstar,xlist,ier,its]

       x0 = x1
       

    xstar = x1
    ier = 1
    return[xstar,xlist,ier,its]


    
x0 = [1,0]
tol = 1e-10
Nmax = 100
xstar,xlist, ier, its = LazyNewton(x0,tol,Nmax)
print(f"The approximate root of f(x) is ", xstar)
print(f'The number of iterations for convergence was {its} iterations')
print(f"Error message: ", ier)