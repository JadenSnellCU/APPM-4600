import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
from numpy.polynomial.laguerre import laggauss

def composite_trapezoidal_rule(a, b, f, N,k):
    delta_x = (b - a) / (N - 1)
    x_values = [a + i*delta_x for i in range(N)]
    function_evaluations = 0
    sum_f_x = 0
    for x in x_values:
        sum_f_x += f(x, k)
        function_evaluations += 1
    result = delta_x * ((f(a, k) + f(b, k)) / 2 + sum_f_x )
    function_evaluations += 2  
    return result, function_evaluations
    

def test_function(x,k):
    return x**(k-1) * np.exp(-x)

def truncation_point(k, tol = 1e-6):
    x = 0.5
    while True:
        if test_function(x,k) < tol:
            return x
        x+=0.1

for k in [2, 4, 6, 8, 10,20]:
    a, b = 0, truncation_point(k) 
    N = int(b)*100  # Number of intervals
    numeric_integral,neval = composite_trapezoidal_rule(a, b, test_function, N,k)
    gamma_value = gamma(k)
    print(f"Gamma({k}) = {gamma_value}, Numeric = {round(numeric_integral,6)}, Relative Error = {abs(gamma_value - numeric_integral)/gamma_value}, Function Evaluations = {neval}")
    print("Truncation point:", round(b,2))
    
#Python quad
print("\nPython Quadrature")
for k in [2, 4, 6, 8, 10, 20]:
    a, b = 0, truncation_point(k)
    integral_result = quad(test_function, a, b, args=(k),full_output=1)
    gamma_value = gamma(k)
    numeric_integral = integral_result[0]
    neval = integral_result[2]['neval']
    print(f"Gamma({k}) = {gamma_value}, Numeric Quadrature Routine = {round(numeric_integral, 6)}, Relative Error = {abs(gamma_value - numeric_integral) / gamma_value}, Function Evaluations = {neval}")
    print("Truncation point:", round(b, 2))


def gamma_approx(t, n):
    # Obtain the n points and weights for Gauss-Laguerre quadrature
    x, w = laggauss(n)
    # Calculate the approximate integral for Gamma(t)
    integral_approx = np.sum(w * x**(t-1))
    return integral_approx

t_values = [2, 3, 4, 5, 6]  # Example values for t


for k in [2, 4, 6, 8, 10, 20]:
    approx = gamma_approx(k, 10)
    exact = gamma(k)
    print(f"Gamma({k}) = {exact}, Approximation = {approx}, Relative Error = {abs(exact - approx)/exact}")