import math


x = 9.999999995000000 * 10**(-10)

#Original
y = math.exp(x)
result_1 = y - 1 

#Approx
y_approx = 1 +  x
result_2 = y_approx - 1

print(f"Original result: {result_1:.16f}")
print(f"Taylor Series result: {result_2:.16f}")