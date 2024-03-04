import numpy as np
import math
import time
from numpy.linalg import inv
from numpy.linalg import norm
import matplotlib.pyplot as plt;

def cent_diff(f, h, x, j):
    # Modify to use index j for adaptive hj
    perturb = np.zeros_like(x)
    perturb[j] = h 
    df = (f(x + perturb) - f(x - perturb)) / (2 * h)
    return df

def evalF(x):
    F = np.zeros(2)
    F[0] = 4*x[0]**2 + x[1]**2 - 4
    F[1] = x[0] + x[1] - np.sin(x[0] - x[1])
    return F

def evalJ(x, h_factors):
    h = [h_factors[i] * abs(x[i]) for i in range(len(x))]
    J = np.array([[cent_diff(lambda x: evalF(x)[0], h[0], x, 0), cent_diff(lambda x: evalF(x)[0], h[1], x, 1)],
                  [cent_diff(lambda x: evalF(x)[1], h[0], x, 0), cent_diff(lambda x: evalF(x)[1], h[1], x, 1)]])
    return J


def LazyNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    xlist = np.zeros((Nmax+1,len(x0)))
    xlist[0] = x0
    h_factors = np.array([1e-3, 1e-3])
    J = evalJ(x0,h_factors)
    
    for its in range(Nmax):
       if its % 3 == 0:
           h_factors /=2
           J = evalJ(x0,h_factors) 
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


x0 = [1.1,0.2]
tol = 1e-10
Nmax = 100
xstar,xlist, ier, its = LazyNewton(x0,tol,Nmax)
print(f"The approximate root of f(x) is ", xstar)
print(f'The number of iterations for convergence was {its} iterations')
print(f"Error message: ", ier)