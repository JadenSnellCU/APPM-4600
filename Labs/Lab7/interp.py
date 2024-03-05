import numpy as np
import matplotlib.pyplot as plt

def dividedDiffTable(x, y, n):
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = (y[j][i - 1] - y[j + 1][i - 1]) / (x[j] - x[i + j])
    return y

def eval_lagrange(x, xint, yint, N):
    lj = np.ones(N + 1)
    for count in range(N + 1):
        for jj in range(N + 1):
            if jj != count:
                lj[count] *= (x - xint[jj]) / (xint[count] - xint[jj])
    yeval = sum(yint[jj] * lj[jj] for jj in range(N + 1))
    return yeval

def evalDDpoly(x, xint, y, N):
    ptmp = np.zeros(N + 1)
    ptmp[0] = 1.
    for j in range(N):
        ptmp[j + 1] = ptmp[j] * (x - xint[j])
    yeval = sum(y[0][j] * ptmp[j] for j in range(N + 1))
    return yeval

def monomial(xeval, xint, yint, N):
    v = np.vander(xint, N + 1, increasing=True)
    a = np.linalg.solve(v, yint)
    yeval = np.polyval(a[::-1], xeval)
    return yeval

def evaluate_and_plot(N_values, function, function_label):
    a, b = -1, 1
    Neval = 1000
    xeval = np.linspace(a, b, Neval + 1)

    for N in N_values:
        xint = np.linspace(a, b, N + 1)
        yint = function(xint)

        y = np.zeros((N + 1, N + 1))
        for j in range(N + 1):
            y[j][0] = yint[j]
        y = dividedDiffTable(xint, y, N + 1)

        yeval_mono = monomial(xeval, xint, yint, N)
        yeval_l = np.array([eval_lagrange(x, xint, yint, N) for x in xeval])
        yeval_dd = np.array([evalDDpoly(x, xint, y, N) for x in xeval])

        fex = function(xeval)

        plt.figure(figsize=(14, 7))
        plt.subplot(1, 2, 1)
        plt.plot(xeval, fex, 'k-', label='Exact')
        plt.plot(xeval, yeval_mono, 'g--', label='Monomial')
        plt.plot(xeval, yeval_l, 'r:', label='Lagrange')
        plt.plot(xeval, yeval_dd, 'b-.', label='Newton DD')
        plt.title(f'Approximations for N={N}')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.semilogy(xeval, np.abs(fex - yeval_mono), 'g--', label='Monomial')
        plt.semilogy(xeval, np.abs(fex - yeval_l), 'r:', label='Lagrange')
        plt.semilogy(xeval, np.abs(fex - yeval_dd), 'b-.', label='Newton DD')
        plt.title(f'Log Absolute Error for N={N}')
        plt.legend()

        plt.suptitle(f'{function_label} Function Interpolation and Error Analysis')
        plt.show()

        if np.max(yeval_mono) >= 100:
            print(f"Stopping: Maximum value of p(x) reached 100 for N={N}")
            break

# Function definitions
f = lambda x: 1 / (1 + (10 * x) ** 2)
g = lambda x: np.sinc(5 * x)

# Evaluation and plotting
N_values_initial = range(2, 11)  # For question 2
N_values_continuation = range(11, 21)  #for question 3
evaluate_and_plot(N_values_initial, f, "f(x) = 1 / (1 + (10 * x) ^ 2)")
evaluate_and_plot(N_values_continuation, f, "Continuation for f(x) = 1 / (1 + (10 * x) ^ 2)")
evaluate_and_plot(N_values_initial, g, "g(x) = sinc(5x)")
