import numpy as np
from scipy.integrate import quad

def composite_trapezoidal_rule(a, b, f, N):
    delta_x = (b - a) / (N - 1)
    x_values = [a + i*delta_x for i in range(N)]
    sum_f_x = sum(f(x) for x in x_values[1:-1])
    return delta_x * ((f(a) + f(b)) / 2 + sum_f_x)

def composite_simpsons_rule(a, b, f, N):
    if N % 2 == 1:  
        N+= 1
    h = (b - a) / N
    x = np.linspace(a, b, N+1)
    y = f(x)
    return h / 3 * (y[0] + 2 * np.sum(y[2:N:2]) + 4 * np.sum(y[1:N:2]) + y[N])

def test_function(x):
    return 1/(1+x**2)

a = -5
b = 5
N_trap = 1291
N_simp = 108
actual = 2 * np.arctan(5)
tol4=1e-4
tol6=1e-6

result_4, err_4,info_4 = quad(test_function, a, b, epsabs=1e-4,full_output=1)
result_6, err_4,info_6 = quad(test_function, a, b, epsabs=1e-6,full_output=1)

iter4 = info_4['neval']
iter6 = info_6['neval']

trapezoidal_result = composite_trapezoidal_rule(a, b, test_function, N_trap)
print("Composite Trapezoidal Rule Result:", round(trapezoidal_result,10))
print(f"Trapezoidal Rule Error for N = {N_trap}:",np.abs(trapezoidal_result - actual) )

simpsons_result = composite_simpsons_rule(a, b, test_function, N_simp)
print("\nComposite Simpson's Rule Result:", round(simpsons_result,10))
print(f"Simpson's Rule Error for N = {N_simp}:",np.abs(simpsons_result - actual) )

for i in range(2,500):
    trapezoidal_result = composite_trapezoidal_rule(a, b, test_function, i)
    if np.abs(trapezoidal_result-actual) < tol4 :
        iterT4=i
        break

for i in range(2,50):
    simpsons_result =composite_simpsons_rule(a, b, test_function, 2*i)
    if np.abs(simpsons_result-actual) < tol4 :
        iterS4=2*i
        break
    
for i in range(2,1000):
    trapezoidal_result = composite_trapezoidal_rule(a, b, test_function, i)
    if np.abs(trapezoidal_result-actual) < tol6 :
        iterT6=i
        break

for i in range(2,50):
    simpsons_result =composite_simpsons_rule(a, b, test_function, 2*i)
    if np.abs(simpsons_result-actual) < tol6 :
        iterS6=2*i
        break
    
print("\nTrapezidal rule function evaluations for tol = 10^-4:", iterT4)
print("Trapezoidal rule function evaluations for tol = 10^-6:", iterT6)

print("\nSimpson's rule function evaluations for tol = 10^-4:", iterS4)
print("Simpson's rule function evaluations for tol = 10^-:6", iterS6)

print("\nQuadratue function evaluations for tol = 10^-4:", iter4)
print("Quadratue function evaluations for tol = 10^-6:", iter6)