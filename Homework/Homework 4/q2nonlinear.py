import numpy as np
from numpy import random as rand
import time
import math
from scipy import io, integrate, linalg, signal
from scipy.linalg import lu_factor, lu_solve
import matplotlib.pyplot as plt
from numpy.linalg import norm

f = lambda x,y: x**2 + y**2 - 4
g = lambda x,y: np.exp(x) + y - 1
Nmax = 100
tol = 1e-10

#Evaluation of the Jacobian
def evalJ(x):
    J = np.array([[2*x[0], 2*x[1]],[np.exp(x[0]), 1]])
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


#Lazy Newton
def LazyNewton(x_0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    xlist = np.zeros((Nmax+1,len(x_0)))
    xlist[0] = x_0
    for its in range(Nmax):
        if its % 3 == 0:
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

#Broyden's Method
def Broyden(x_0, tol, Nmax):
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, xlist = list of approximations, ier = error message, its = num its'''
    xlist = np.zeros((Nmax+1, len(x_0)))
    xlist[0] = x_0
    B = evalJ(x_0)  # Initial approximation to the Jacobian

    for its in range(Nmax):
        F = np.array([f(x_0[0], x_0[1]), g(x_0[0], x_0[1])])
        s = np.linalg.solve(B, -F)
        x1 = x_0 + s
        xlist[its+1] = x1

        if norm(s) < tol*norm(x_0):
            xstar = x1
            ier = 0
            return [xstar, xlist, ier, its]

        # Update the approximation to the Jacobian
        y = np.array([f(x1[0], x1[1]), g(x1[0], x1[1])]) - F
        z = y - B @ s
        B += np.outer(z, s) / (s @ s)

        x_0 = x1

    xstar = x1
    ier = 1
    return [xstar, xlist, ier, its]

#Broyden's Method
def Broyden_id(x_0, tol, Nmax):
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, xlist = list of approximations, ier = error message, its = num its'''
    xlist = np.zeros((Nmax+1, len(x_0)))
    xlist[0] = x_0
    B = np.array([[1.0,0.0],[0.0,1.0]]) # Use Identity matrix

    for its in range(Nmax):
        F = np.array([f(x_0[0], x_0[1]), g(x_0[0], x_0[1])])
        s = np.linalg.solve(B, -F)
        x1 = x_0 + s
        xlist[its+1] = x1

        if norm(s) < tol*norm(x_0):
            xstar = x1
            ier = 0
            return [xstar, xlist, ier, its]

        # Update the approximation to the Jacobian
        y = np.array([f(x1[0], x1[1]), g(x1[0], x1[1])]) - F
        z = y - B @ s
        B += np.outer(z, s) / (s @ s)

        x_0 = x1

    xstar = x1
    ier = 1
    return [xstar, xlist, ier, its]

#i
x_0 = [1,1]
print("\nInitial guess of x = 1, y = 1\n")

# Call Newton's method 
xstar_newton, xlist_newton, ier_newton, its_newton = Newton(x_0, tol, Nmax)

# Output the results for Newton's method
if ier_newton == 0:
    print(f"Newton's method converged after {its_newton} iterations: x = {xstar_newton[0]}, y = {xstar_newton[1]}")
else:
    print("Newton's method did not converge")


# Call the Lazy Newton's method function
xstar_lazy_newton, xlist_lazy_newton, ier_lazy_newton, its_lazy_newton = LazyNewton(x_0, tol, Nmax)

# Output the results for Lazy Newton's method
if ier_lazy_newton == 0:
    print(f"Lazy Newton's method converged after {its_lazy_newton} iterations: x = {xstar_lazy_newton[0]}, y = {xstar_lazy_newton[1]}")
else:
    print("Lazy Newton's method did not converge")


# Call the Broyden's method function
xstar_broyden, xlist_broyden, ier_broyden, its_broyden = Broyden(x_0, tol, Nmax)

# Output the results for Broyden's method
if ier_broyden == 0:
    print(f"Broyden's method converged after {its_broyden} iterations: x = {xstar_broyden[0]}, y = {xstar_broyden[1]}")
else:
    print("Broyden's method did not converge")




#ii
x_1 = [1,-1]
print("\nInitial guess of x = 1, y = -1\n")
# Call Newton's method 
xstar_newton, xlist_newton, ier_newton, its_newton = Newton(x_1, tol, Nmax)

# Output the results for Newton's method
if ier_newton == 0:
    print(f"Newton's method converged after {its_newton} iterations: x = {xstar_newton[0]}, y = {xstar_newton[1]}")
else:
    print("Newton's method did not converge")


# Call the Lazy Newton's method function
xstar_lazy_newton, xlist_lazy_newton, ier_lazy_newton, its_lazy_newton = LazyNewton(x_1, tol, Nmax)

# Output the results for Lazy Newton's method
if ier_lazy_newton == 0:
    print(f"Lazy Newton's method converged after {its_lazy_newton} iterations: x = {xstar_lazy_newton[0]}, y = {xstar_lazy_newton[1]}")
else:
    print("Lazy Newton's method did not converge")


# Call the Broyden's method function
xstar_broyden, xlist_broyden, ier_broyden, its_broyden = Broyden(x_1, tol, Nmax)

# Output the results for Broyden's method
if ier_broyden == 0:
    print(f"Broyden's method converged after {its_broyden} iterations: x = {xstar_broyden[0]}, y = {xstar_broyden[1]}")
else:
    print("Broyden's method did not converge")
    
    

#iii
x_2 = [0,0]
print("\nInitial guess of x = 0, y = 0\n")
'''
# Call Newton's method 
xstar_newton, xlist_newton, ier_newton, its_newton = Newton(x_2, tol, Nmax)

# Output the results for Newton's method
if ier_newton == 0:
    print(f" Newton's method converged after {its_newton} iterations: x = {xstar_newton[0]}, y = {xstar_newton[1]}")
else:
    print("Newton's method did not converge")

'''
'''
# Call the Lazy Newton's method function
xstar_lazy_newton, xlist_lazy_newton, ier_lazy_newton, its_lazy_newton = LazyNewton(x_2, tol, Nmax)

# Output the results for Lazy Newton's method
if ier_lazy_newton == 0:
    print(f"Lazy Newton's method converged after {its_lazy_newton} iterations: x = {xstar_lazy_newton[0]}, y = {xstar_lazy_newton[1]}")
else:
    print("Lazy Newton's method did not converge")

'''
# Call the Broyden's method function (use identity matrix because of singular jacobian) 
xstar_broyden, xlist_broyden, ier_broyden, its_broyden = Broyden_id(x_2, tol, Nmax)

# Output the results for Broyden's method
if ier_broyden == 0:
    print(f"Broyden's method converged after {its_broyden} iterations: x = {xstar_broyden[0]}, y = {xstar_broyden[1]}")
else:
    print("Broyden's method did not converge")