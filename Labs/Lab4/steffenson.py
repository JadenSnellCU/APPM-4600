import numpy as np

# test functions
f1 = lambda x: (10/(x+4))**0.5 # fixed point is f1 = 1.3652300134140976

Nmax = 100
tol = 1e-11

# refactored subroutine
def fixedpt(f, x0, tol, Nmax):
    x = np.zeros(Nmax+1)
    x[0] = x0
    ier = 1
    for count in range(1, Nmax+1):
        x[count] = f(x[count-1])
        if abs(x[count] - x[count-1]) < tol:
            ier = 0
            break
    return x[:count+1], ier,count

#Calculate a and lambda 
def estimate_order(x, tol):
    m = len(x) - 1  
    if m < 2:  
        return None, None
    
    # Using the formula for order of convergence and asymptotic error constant
    a = np.log(abs(x[m] - x[m-1]) / abs(x[m-1] - x[m-2])) / np.log(abs(x[m-1] - x[m-2]) / abs(x[m-2] - x[m-3]))
    err_lambda = abs(x[m] - x[m-1]) / abs(x[m-1] - x[m-2])**a
    return a, err_lambda

#The needed inputs for steffenson's method is the sequence from fixed point  and it will output a sequence used with the approximations 
def steffensons(g, p0, tol, Nmax):
   count = 0
   p1 = -1
   while count < Nmax:
       count += 1
       a = p0
       b = g(a)
       c = g(b)
       p1 = a - ((b - a) ** 2) / (c - 2 * b + a)
       if abs(p1 - p0) < tol:
           return [p1, 0]
       else:
           p0 = p1

   return [p1, 1]
