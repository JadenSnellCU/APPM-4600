import numpy as np
from scipy.special import gamma

def composite_trapezoidal_rule(a, b, f, N,k):
    delta_x = (b - a) / (N - 1)
    x_values = [a + i*delta_x for i in range(N)]
    sum_f_x = sum(f(x,k) for x in x_values[1:-1])
    return delta_x * ((f(a,k) + f(b,k)) / 2 + sum_f_x)

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
    numeric_integral = composite_trapezoidal_rule(a, b, test_function, N,k)
    gamma_value = gamma(k)
    print(f"Gamma({k}) = {gamma_value}, Numeric = {round(numeric_integral,6)}, Relative Error = {abs(gamma_value - numeric_integral)/gamma_value}")
    print("Truncation point:", round(b,2))