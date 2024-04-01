import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def driver():
    f = lambda x: np.sin(9*x)
    a, b = 0, (2*np.pi)/9

    num = [5, 10, 20, 40]
    for Nint in num:
        xint = np.linspace(a, b, Nint)
        yint = f(xint)

        # Use the scipy CubicSpline to create a periodic spline
        cs = CubicSpline(xint, yint, bc_type='periodic')

        xeval = np.linspace(a, b, 1000)
        yeval = cs(xeval)

        # Plot the function and the spline
        plt.figure()
        plt.plot(xeval, f(xeval), label='True function')
        plt.plot(xeval, yeval, label=f'Spline with {Nint} points')
        plt.legend()
        plt.show()

        # Error analysis
        error = np.log10(np.abs(f(xeval) - yeval))
        plt.figure()
        plt.plot(xeval, error, label=f'Log Error with {Nint} points')
        plt.legend()
        plt.show()

driver()
