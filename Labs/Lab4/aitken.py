import numpy as np

#p_{n+1}^2 -2p_{n+1}*p = p_{n+2}*p_n - p*p_{n+2} - p*p_n

#p(p_{n+2}+p_n-2p_{n+1}) = p_{n+2}*p_n - p_{n+1}^2

# p =  ( p_{n+2}*p_n - p_{n+1}^2)/(p_{n+2}+p_n-2p_{n+1})


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

def aitken(approx, tol, Nmax):
    x = np.zeros(len(approx))
    ier = 1
    if len(approx)< 3:
        return None, None
    for count in range(0,len(approx)-1):
        x[count] = approx[count]-((approx[count+1]-approx[count])**2)/(approx[count+2] - 2*approx[count+1] + approx[count])
        if abs(x[count] - x[count-1]) < tol:
            ier = 0
            break
    return x[:count+1], ier, count
# test f1
x0 = 1.5
x_approximations, ier, count = fixedpt(f1, x0, tol, Nmax)
a1, err_lambda1 = estimate_order(x_approximations, tol)
print('The sequence for f1 is:', x_approximations)
print('Order of convergence for f1:', round(a1,10))
print('Asymptotic error constant for f1:', err_lambda1)
print('Number of iterations: ', count)
print('Error message reads:', ier)
print('\n')

accel_approx, ier, count = aitken(x_approximations,tol,Nmax)
a2, err_lambda2 = estimate_order(accel_approx, tol)
print('The sequence for acceleration method is:', accel_approx)
print('Order of convergence for acceleration method:', round(a2,10))
print('Asymptotic error constant for acceleration method:', err_lambda2)
print('Number of iterations: ', count)
print('Error message reads:', ier)