# import libraries
import numpy as np


# define routines
def bisection(f,a,b,tol,Nmax):
    '''
    Inputs:
      f,a,b       - function and endpoints of initial interval
      tol, Nmax   - bisection stops when interval length < tol
                  - or if Nmax iterations have occured
    Returns:
      astar - approximation of root
      ier   - error message
            - ier = 1 => cannot tell if there is a root in the interval
            - ier = 0 == success
            - ier = 2 => ran out of iterations
            - ier = 3 => other error ==== You can explain
    '''

    '''     first verify there is a root we can find in the interval '''
    fa = f(a); fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

    ''' verify end point is not a root '''
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    while (count < Nmax):
      c = 0.5*(a+b)
      fc = f(c)

      if (fc ==0):
        astar = c
        ier = 0
        return [astar, ier]

      if (fa*fc<0):
         b = c
      elif (fb*fc<0):
        a = c
        fa = fc
      else:
        astar = c
        ier = 3
        return [astar, ier]

      if (abs(b-a)<tol):
        astar = a
        ier =0
        return [astar, ier]
      
      count = count +1

    astar = a
    ier = 2
    return [astar,ier] 

#Exercises
f = lambda x: x**3-x**2
Nmax = 100
tol = 10e-5

#4.1  
#a 
print("Exercise 4.1(a)")
a = 0.5
b = 2


[astar,ier] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',astar)
print(f'the error message reads  for the interval [{a},{b}] is',ier)
print("\n")

#b
print("Exercise 4.1(b)")
a = -1
b = 0.5


[astar,ier] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',astar)
print(f'the error message reads  for the interval [{a},{b}] is',ier)
print('As seen by the results above, bisection cannot find a root for [-1,0.5] as both values are negative for f(x)')
print("\n")


#c  
print("Exercise 4.1(c)")
a = -1
b = 2

[astar,ier] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',astar)
print(f'the error message reads  for the interval [{a},{b}] is',ier)
print("\n")


#4.2
#a
print("Exercise 4.2(a)")
f = lambda x: (x-1)*(x-3)*(x-5)
a = 0
b = 2.4
[astar,ier] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',astar)
print(f'the error message reads  for the interval [{a},{b}] is',ier)
print("\n")

#b
print("Exercise 4.2(b)")
f = lambda x: (x-1)**2 * (x-3) 
a=0 
b=2
[astar,ier] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',astar)
print(f'the error message reads  for the interval [{a},{b}] is',ier)
print("\n")

#c
print("Exercise 4.2(c) for a = 0 and b = 0.1")
f = lambda x: np.sin(x)
a = 0
b = 0.1
[astar,ier] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',astar)
print(f'the error message reads  for the interval [{a},{b}] is',ier)
print("Due to the both f(a) and f(b) being positive, it's impossible to find the root through bisection.")
print("\n")


f = lambda x: np.sin(x)
a = 0.5
b = 3* np.pi / 4
print(f"Exercise 4.2(c) for a = 0.5 and b = {b} ")
[astar,ier] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',astar)
print(f'the error message reads  for the interval [{a},{b}] is',ier)
print("Due to the both f(a) and f(b) being positive, it's impossible to find the root through bisection.")
print("\n")