import numpy as np
import matplotlib.pyplot as plt
import random

#Problem a
# Creating the vector t
t = np.linspace(0, np.pi, 30)

# Creating the vector y as the cosine of t
y = np.cos(t)

# Calculating the sum S
S = np.sum(t * y)

# Printing the result
print("The sum is " + str(S))


#Problem 2

#Initializing parameters
R = 1.2
delta_r = 0.1
f = 15
p = 0
theta = np.linspace(0,2*np.pi, 100)

#Initialize equations for first figure
x = R*(1+delta_r*np.sin(f*theta +p))*np.cos(theta)
y = R*(1+delta_r*np.sin(f*theta +p))*np.sin(theta)

#Plot
fig1,ax1 = plt.subplots()
ax1.plot(x,y)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
plt.title("Parametric Curve")
plt.show()

#Reinitialize parameters
num_curves = 10
delta_r = 0.05

#Create the curves
fig2, ax2 = plt.subplots()
for i in range(1, num_curves + 1):
    R = i
    f = 2 + i

    # Generate a random value for p
    p = random.uniform(0, 2)

    x = R*(1 + delta_r*np.sin(f*theta + p))*np.cos(theta)
    y = R*(1 + delta_r*np.sin(f*theta + p))*np.sin(theta)

    ax2.plot(x, y, label=f'Curve {i}')

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.legend()
plt.title("Wavy Bands")
plt.show()






