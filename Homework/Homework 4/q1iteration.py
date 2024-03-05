import numpy as np
from numpy import random as rand
import time
import math
from scipy import io, integrate, linalg, signal
from scipy.linalg import lu_factor, lu_solve
import matplotlib.pyplot as plt
from numpy.linalg import norm



f = lambda x,y: 3*x**2 - y**2
g = lambda x,y: 3*x*(y**2)-x**3-1

coeff_mat = np.array([[1/6,1/18],
                     [0,1/6]])
n_max = 100

tol = 1e-10

def iterate_scheme(f,g,coeff_mat,n_max,tol):
    """Using the iteration scheme discussed in Problem 1
    err_message meaning: 
    err_msg = 0 : Root was found
    err_msg = 1: Error, did not converge
    
    """
    x_0 = 1 ; y_0 = 1
    err_msg = 0
    for i in range(n_max):
        F = f(x_0,y_0)
        G = g(x_0,y_0)
        
        x_update, y_update = coeff_mat @ np.array([F,G])
        
        x_n1 , y_n1 = x_0 - x_update, y_0 - y_update
        
        if np.abs(x_n1 - x_0) < tol and np.abs(y_n1 - y_0) < tol:
            return x_n1, y_n1, err_msg, i+1
        
        x_0,y_0 = x_n1, y_n1
    err_msg = 1
    return x_n1, y_n1, err_msg, n_max
    
x_sol, y_sol, ier, count = iterate_scheme(f,g,coeff_mat,n_max,tol)
if ier == 0:
    print(f"This iterative scheme converged after {count} iterations: x = {x_sol}, y = {y_sol}")
else:
    print("Did not converge")
     

#Evaluation of the Jacobian
def evalJ(x):
    J = np.array([[6*x[0], -2*x[1]],[3*x[1]**2-3*x[0]**2, 6*x[0]*x[1]]])
    return J


#Newtons Method
def Newton(x_0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    xlist = np.zeros((Nmax+1,len(x_0)))
    xlist[0] = x_0
    for its in range(Nmax):
       J = evalJ(x_0)
       F = [f(x_0[0],x_0[1]),g(x_0[0],x_0[1])]

       x1 = x_0 - np.linalg.solve(J,F)
       xlist[its+1]=x1

       if (norm(x1-x_0) < tol*norm(x_0)):
           xstar = x1
           ier =0
           return[xstar, xlist,ier, its]

       x_0 = x1

    xstar = x1
    ier = 1
    return[xstar,xlist,ier,its]

x_0 = [1,1]

sols, sol_list, ier, count = Newton(x_0,tol,n_max)
if ier == 0:
    print(f"Newton's method converged after {count} iterations: x = {sols[0]}, y = {sols[1]}")
else:
    print("Did not converge")