import numpy as np
import matplotlib.pyplot as plt
#1-3
x = np.linspace(0,10,100)
y=np.arange(0,10,0.1)
print('The first three entries of x are', x[:3])

#4-6
w = 10**(-np.linspace(1,10,10))
s = 3*w
x = np.linspace(1,len(w), len(w))

plt.semilogy(x,w, label = 'W vector')
plt.semilogy(x,s, label = 'S vector')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()


