import numpy as np


def forw_diff(f,h,x0):
    for h_value in h:
        df = (f(x0+h_value)-f(x0))/h_value
    return df

def cent_diff(f,h,x0):
    for h_value in h:
        df = (f(x0+h_value)-f(x0-h_value))/(2*h_value)
    return df

def estimate_order(x):
    m = len(x) - 1  
    if m < 2:  
        return None, None
    
    # Using the formula for order of convergence and asymptotic error constant
    a = np.log(abs(x[m] - x[m-1]) / abs(x[m-1] - x[m-2])) / np.log(abs(x[m-1] - x[m-2]) / abs(x[m-2] - x[m-3]))
    err_lambda = abs(x[m] - x[m-1]) / abs(x[m-1] - x[m-2])**a
    return a, err_lambda

#Function
f = lambda x: np.cos(x)

#Step
h = [0.01 / (2**(np.arange(0,10)))]



#Point of Evaluation
x0 = np.pi/2


#Forward Difference
df = forw_diff(f,h,x0)
a, err_lambda = estimate_order(df)
print(f"The derivative at point x = {x0} approximated through forward difference is " , df[1])
print(f"The order of convergence of forward difference is ",  round(a,2))
print('\n')

#Centered Difference
df = cent_diff(f,h,x0)
a, err_lambda = estimate_order(df)
print(f"The derivative at point x = {x0} approximated through centered difference is " , df[-1])
print(f"The order of convergence of centered difference is ", round(a,2))
print('\n')
