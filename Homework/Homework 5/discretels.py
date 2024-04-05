import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 1, 2, 3])
y = np.array([1, 4, 2, 6])

# Constructing the G matrix
G = np.vstack((np.ones_like(x), x)).T

# Constructing the b vector
b = y

GT = G.T
GTG = np.dot(GT, G)
GTb = np.dot(GT, b)


a = np.linalg.solve(GTG, GTb)
print(a)  



weights = np.array([1, 4, 9, 6])

# Constructing the D matrix
D = np.diag(np.sqrt(weights))

# Calculating the weighted normal equation components
MtD = G.T @ D
MtDM = MtD @ G
MtDy = MtD @ y



a_weighted = np.linalg.solve(MtDM, MtDy)
print(a_weighted) 


x_plot = np.linspace(0, 3, 100)

# Compute the y values of the fit lines
y_fit_unweighted = a[0] + a[1] * x_plot
y_fit_weighted = a_weighted[0] + a_weighted[1] * x_plot

plt.scatter(x, y, color='purple', label='Data points')
plt.plot(x_plot, y_fit_unweighted, label='Unweighted least squares', color='blue')
plt.plot(x_plot, y_fit_weighted, label='Weighted least squares', color='red')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Comparison of Least Squares Fits')
plt.legend()

plt.show()