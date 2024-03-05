import numpy as np
from numpy import random as rand
import time
import math
from scipy import io, integrate, linalg, signal
from scipy.linalg import lu_factor, lu_solve
import matplotlib.pyplot as plt
from numpy.linalg import norm


f = lambda x,y,z: x + np.cos(x*y*z) - 1 

g = lambda x,y,z : (1-x)**(1/4) + y + 0.05*z**2 - 0.15*z - 1

h = lambda x,y,z: -x**2 - 0.1*y **2 + 0.01*y + z - 1

Nmax = 100
tol = 1e-6

#Evaluation of the Jacobian
def evalJ(x):
    '''
    inputs: x = nonlinear system of equations
    outputs: J = jacobian of system
    '''
    J = np.array([[1 - x[1] * x[2] * np.sin(x[0]*x[1]*x[2]), - x[0]*np.sin(x[0]*x[1]*x[2]), - x[2]*np.sin(x[0]*x[1]*x[2])],
                  [-(1/4)/((1-x[0])**(3/4)), 1, .1*x[2] - 0.15],
                  [-2*x[0], -0.2 * x[1] + 0.01, 1]])
    return J

#Newtons Method
def Newton(x_0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx sol to system, ier = error message, its = num its'''

    for its in range(Nmax):
        J = evalJ(x_0)
        F = [f(x_0[0],x_0[1],x_0[2]),g(x_0[0],x_0[1], x_0[2]), h(x_0[0],x_0[1],x_0[2])]

        x1 = x_0 - np.linalg.solve(J,F)

        if (norm(x1-x_0) < tol*norm(x_0)):
           xstar = x1
           ier =0
           return xstar,ier, its

        x_0 = x1
    xstar = x1
    ier = 1
    return xstar,ier,its

#Steepest Descent
def steepest_descent(x_0,tol, Nmax):
    x = x_0
    for its in range(Nmax):
        J = evalJ(x)  
        F = np.array([f(x[0], x[1], x[2]), g(x[0], x[1], x[2]), h(x[0], x[1], x[2])])
        grad_F = J.T @ F #Find the gradient

        # Determine the descent direction for line search func
        d = -grad_F

        # Define the objective function for line search
        res = lambda x: np.sum([f(*x)**2, g(*x)**2, h(*x)**2])
        grad = lambda x: 2 * J.T @ np.array([f(*x), g(*x), h(*x)])

        # Find alpha
        alpha = line_search(res, grad, x, d)

        
        x_next = x + alpha * d

       
        if norm(x_next - x) < tol:
            ier = 0
            return x_next, ier, its 

        x = x_next  

    return x,ier, Nmax  

def line_search(f, grad_f, x, dir, a=1, rho=0.5, c=1e-4):
    alpha = a
    while f(x+ alpha *dir) > f(x) + c * alpha * np.dot(grad_f(x),dir):
        alpha*=rho
        if alpha < 1e-10:
            break
    return alpha
        
x_0 = [0,1,1]
print("\nInitial guess of x = 0, y = 1, z = 1\n")

# Call Newton's method 
xstar_newton, ier_newton, its_newton = Newton(x_0, tol, Nmax)

# Output the results for Newton's method
if ier_newton == 0:
    print(f"\nNewton's method converged after {its_newton} iterations: x = {xstar_newton[0]}, y = {xstar_newton[1]}, z = {xstar_newton[2]}\n")
else:
    print("Newton's method did not converge")
    

# Call Steepest Descent method 
xstar_steep, ier_steep, its_steep = steepest_descent(x_0, tol, Nmax)

# Output the results for Steepest Descent method
if ier_steep == 0:
    print(f"\nSteepest Descent method converged after {its_steep} iterations: x = {xstar_steep[0]}, y = {xstar_steep[1]}, z = {xstar_steep[2]}\n")
else:
    print("Steepest Descent method did not converge")
    

#Hybrid method

hybrid_tol = 5e-2
newt_init, ier, iter = steepest_descent(x_0,hybrid_tol,Nmax)
xstar_newton, ier_newton, its_newton = Newton(newt_init, tol, Nmax)
   
# Output the results for Hybrid method
if ier_newton == 0:
    print(f"\nHybrid method converged after {its_newton + iter} iterations: x = {xstar_newton[0]}, y = {xstar_newton[1]}, z = {xstar_newton[2]}\n")
else:
    print("Hybrid method did not converge")
