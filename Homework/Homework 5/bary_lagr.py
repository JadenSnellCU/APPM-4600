import numpy as np
import matplotlib.pyplot as plt

#Barycentric Lagrange Interpolation Function
def bary_lagr(x,f,z):
    n = len(x)
    w = np.zeros(n)
    for j in range(0,n):
        w[j] = 1
        for i in range(0,n):
            if i != j:
                w[j] = w[j]/(x[j]-x[i])
            else:
                continue
    p = np.zeros(len(z))
    
    for k in range(0, len(z)):
        num = 0
        den = 0
        
        for j in range(0,n):
            if z[k] == x[j]:
                p[k] = f[j]
                break
        if p[k] == 0:
            for j in range(0,n):
                t = w[j]/(z[k]-x[j])
                num = num + t*f[j]
                den = den + t
            p[k] = num/den
    return p


f = lambda x: 1/(1+(16*x)**2)

n_values = [2, 3, 4, 5, 10, 15, 20]  
x_range = np.linspace(-1, 1, 1001) 

#Equispaced
for n in n_values:
    h = 2 / n
    x = np.array([-1 + (i -1)* h for i in range(1,n+2)])
    f_x = f(x)
    
    
    # Interpolate
    p_x = bary_lagr(x, f_x, x_range)
    
    # Plot
    plt.figure()
    plt.plot(x, f_x, 'o', label='Data points')
    plt.plot(x_range, p_x, label='Interpolated polynomial')
    plt.title(f'Barycentric Lagrange Interpolation with n={n} using Equispaced points')
    plt.legend()
    plt.show()

    
    if np.max(p_x) >= 100 or np.min(p_x) <= -100:
        print(f"The Runge phenomenon gets extreme when n={n}")
        break
    
#Chebyshev
for n in n_values:
    x = np.array([np.cos(((2*j+1)*np.pi)/(2*(n+1))) for j in range(n+1)])
    f_x = f(x)
    
    
    # Interpolate
    p_x = bary_lagr(x, f_x, x_range)
    
    # Plot
    plt.figure()
    plt.plot(x, f_x, 'o', label='Data points')
    plt.plot(x_range, p_x, label='Interpolated polynomial')
    plt.title(f'Barycentric Lagrange Interpolation with n={n} using Chebyshev points')
    plt.legend()
    plt.show()


def psi_n(x, nodes):
    """Compute psi_n(x) given a set of nodes."""
    psi = np.ones_like(x)
    for xi in nodes:
        psi *= (x - xi)
    return psi

# For the largest n
n = 20

# Equispaced nodes
h = 2 / n
equispaced_nodes = np.array([-1 + i * h for i in range(n+1)])

# Chebyshev nodes
chebyshev_nodes = np.array([np.cos(((2*j+1)*np.pi)/(2*(n+1))) for j in range(n+1)])

# Compute psi_n(x) for both sets of nodes
psi_equispaced = psi_n(x_range, equispaced_nodes)
psi_chebyshev = psi_n(x_range, chebyshev_nodes)

# Plot log10 of the absolute values
plt.figure(figsize=(10, 6))
plt.plot(x_range, np.log10(np.abs(psi_equispaced)), label='Equispaced Nodes')
plt.plot(x_range, np.log10(np.abs(psi_chebyshev)), label='Chebyshev Nodes')
plt.title('Log10 of |psi_n(x)| for Equispaced and Chebyshev Nodes')
plt.xlabel('x')
plt.ylabel('log10(|psi_n(x)|)')
plt.legend()
plt.show()
