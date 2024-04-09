import numpy as np

def composite_trapezoidal_rule(a, b, f, N):
    delta_x = (b - a) / (N - 1)
    x_values = [a + i*delta_x for i in range(N)]
    sum_f_x = sum(f(x) for x in x_values[1:-1])
    return delta_x * ((f(a) + f(b)) / 2 + sum_f_x)

def composite_simpsons_rule(a, b, f, N):
    if N % 2 == 1:
        raise ValueError("N must be even for Simpson's rule")
    
    delta_x = (b - a) / (N - 1)
    sum_even = sum(f(a + 2*i*delta_x) for i in range(1, (N//2)))
    sum_odd = sum(f(a + (2*i-1)*delta_x) for i in range(1, N//2 + 1))
    return (delta_x / 3) * (f(a) + 2*sum_even + 4*sum_odd + f(b))

def test_function(x):
    return np.sin(1/x)

a = 0.1
b = 2
N = 500

trapezoidal_result = composite_trapezoidal_rule(a, b, test_function, N)
print("Composite Trapezoidal Rule Result:", round(trapezoidal_result,5))

simpsons_result = composite_simpsons_rule(a, b, test_function, N)
print("Composite Simpson's Rule Result:", round(simpsons_result,5))