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
    count = 0
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

   
    while (count < Nmax):
      c = 0.5*(a+b)
      fc = f(c)

      if (fc ==0):
        astar = c
        ier = 0
        return [astar,ier, count] 

      if (fa*fc<0):
         b = c
      elif (fb*fc<0):
        a = c
        fa = fc
      else:
        astar = c
        ier = 3
        return [astar,ier, count] 

      if (abs(b-a)<tol):
        astar = a
        ier =0
        return [astar,ier, count] 
      
      count = count +1

    astar = a
    ier = 2
    return [astar,ier, count] 

# problem 1c 
f = lambda x: 2*x - 1 - np.sin(x)
a = 0.5
b = 1

Nmax = 100
tol = 1e-9

[astar,ier, count] = bisection(f,a,b,tol,Nmax)
print("\nFor Problem 1c")
print('the approximate root is',astar)
print('the error message reads:',ier)
print('the number of iterations is', count)
print('\n')


#problem 2a
f = lambda x: (x-5) ** 9
a = 4.82
b = 5.2
tol = 1e-4
[astar,ier, count] = bisection(f,a,b,tol,Nmax)
print("For Problem 2a")
print('the approximate root is',astar)
print('the error message reads:',ier)
print('the number of iterations is', count)
print('\n')

#problem 2b
f = lambda x: x**9 - 45*x**8 + 900*x**7 - 10500*x**6 + 78750*x**5 - 393750*x**4 + 1312500*x**3 - 2812500*x**2 + 3515625*x**1 - 1953125
a = 4.82
b = 5.2
tol = 1e-4
[astar,ier, count] = bisection(f,a,b,tol,Nmax)
print("For Problem 2b")
print('the approximate root is',astar)
print('the error message reads:',ier)
print('the number of iterations is', count)
print('\n')

#problem 3b
f = lambda x: x**3 + x -4
a = 1
b = 4
tol = 1e-3
[astar,ier, count] = bisection(f,a,b,tol,Nmax)
print("For Problem 3b")
print('the approximate root is',astar)
print('the error message reads:',ier)
print('the number of iterations is', count)
print('\n') 