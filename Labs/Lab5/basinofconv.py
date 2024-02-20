import numpy as np

def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  count = 0
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it,count]
      p0 = p1
      count += 1
  pstar = p1
  info = 1
  return [p,pstar,info,it,count]
        

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
        return [astar, ier,count]
      
      count = count +1

    astar = a
    ier = 2
    return [astar,ier] 


# define routines
def hybrid(f,f_prime, f_dprime,a,b,tol,Nmax):
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

      if (abs(f(c)*f_dprime(c)/(f_prime(c)**2) < 1)):
        x0 = c
        ier =0
        break
    
      count = count +1 
    x = np.zeros(Nmax+1)
    x[0] = x0
    for it in range(Nmax):
        x1 = x0-(f(x0)/f_prime(x0))
        x[it+1] = x1
        if (abs(x1-x0) < tol):
            astar = x1
            ier = 0
            return [astar,ier]
        x0 = x1
    astar = x1
    ier = 1
    
    return [astar,ier] 

f = lambda x: np.exp(x**2+7*x-30)-1
fp = lambda x: (2*x+7)*np.exp(x**2+7*x-30)
fdp = lambda x: (4*(x**2)+28*x+51)*np.exp(x**2+7*x-30)
Nmax = 100
tol = 10e-8
a=2
b=4.5

print("Bisection Method")
[bi_astar,bi_ier,count] = bisection(f,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',bi_astar)
print(f'the error message reads  for the interval [{a},{b}] is',bi_ier)
print(f'the # of iterations [{a},{b}] is',count)
print("\n")

print("Newton Method")
[p,newton_astar,newton_ier,it,count] = newton(f,fp,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',newton_astar)
print(f'the error message reads  for the interval [{a},{b}] is',newton_ier)
print(f'the # of iterations [{a},{b}] is',count)
print("\n")

print("Hybrid Method")
[hybrid_astar,hybrid_ier] = hybrid(f,fp,fdp,a,b,tol,Nmax)
print(f'the approximate root is for the interval [{a},{b}] is',hybrid_astar)
print(f'the error message reads  for the interval [{a},{b}] is',hybrid_ier)
print("\n")