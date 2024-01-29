#Driver Code + Exercise 1
import numpy as np
import numpy.linalg as la
import math
import time
def driver():
    n = 100
    x = np.linspace(0,2*np.pi,n)
    
    f = lambda x: np.cos(x)
    g = lambda x: np.sin(x)
    
    y = f(x)
    w = g(x)  
    
    dp = dotProduct(y,w,n)
    print('the dot product is : ', dp)
    return
def dotProduct(x,y,n):
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    return round(dp,10)
driver()

#Making a matrix/Exercise 
def mat_mult(mat,vec):
    num_rows = len(mat)
    product = np.array(np.zeros(num_rows))
    for i in range(num_rows):
        product[i] = dotProduct(mat[i],vec,len(mat))    
    return product

mat = np.array([[1,2],[3,4]])
vec = np.array([5,6])

st1 = time.time()
for i in range(10000):
    product = mat_mult(mat,vec)
end1 = time.time()

print('The resulting vector is', product)

st2 = time.time()
for i in range(10000):
    np.matmul(mat,vec)
end2 = time.time()
print('Using matmul, the resulting vector is',np.matmul(mat,vec))

time1 = end1 - st1
time2 = end2 - st2
print('Time to run my code 10000 times: ', time1)
print('Time to run the lib code 10000 times', time2)

