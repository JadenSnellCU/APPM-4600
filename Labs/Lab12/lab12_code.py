import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import time

def driver():

     ''' create  matrix for testing different ways of solving a square 
     linear system'''

     '''' N = size of system'''
     N_values = [100, 500, 1000, 2000, 4000, 5000]
     for N in N_values:
          ''' Right hand side'''
          b = np.random.rand(N,1)
          A = np.random.rand(N,N)
     
          start_time_lu_1 = time.time()
          lu,piv = scila.lu_factor(A)
          end_time_lu_1 = time.time()
          
          start_time_lu_2 = time.time()
          x = scila.lu_solve((lu,piv),b)
          end_time_lu_2 = time.time()
          
          test = np.matmul(A,x)
          residual = la.norm(test-b)
          
          print(f"Residual using LU Factorization for N = {N}:", residual)
          print("Solution Time - Part 1", end_time_lu_1-start_time_lu_1)
          print("Solution Time - Part 2", end_time_lu_2-start_time_lu_2)
          print("Solution Time - Total", end_time_lu_2-start_time_lu_1)
          print("\n")
          
          start_time = time.time()
          x_direct = scila.solve(A, b)
          end_time = time.time()
          
          test_direct = np.matmul(A, x_direct)
          residual_direct = la.norm(test_direct - b)
          print(f"Residual using direct solve for N = {N}:", residual_direct)
          print("Solution Time", end_time-start_time)
          print('\n')
          
     
     ''' Create an ill-conditioned rectangular matrix '''
     N = 10
     M = 5
     A = create_rect(N,M)     
     b = np.random.rand(N,1)
     
     def solve_normal_equation(A, b):
          A_transpose = np.transpose(A)
          x = np.linalg.solve(np.dot(A_transpose, A), np.dot(A_transpose, b))
          return x

     def solve_qr_factorization(A, b):
          Q, R = np.linalg.qr(A)
          x = np.linalg.solve(R, np.dot(Q.T, b))
          return x

     A = create_rect(N, M)
     b = np.random.rand(N)

     x_normal_equation = solve_normal_equation(A, b)
     x_qr_factorization = solve_qr_factorization(A, b)

     residual_normal_equation = la.norm(np.dot(A, x_normal_equation) - b)
     residual_qr_factorization = la.norm(np.dot(A, x_qr_factorization) - b)

     print("Residual norm using normal equation:", residual_normal_equation)
     print("Residual norm using QR factorization:", residual_qr_factorization)
          


     
def create_rect(N,M):
     ''' this subroutine creates an ill-conditioned rectangular matrix'''
     a = np.linspace(1,15,M)
     d = 10**(-a)
     
     D2 = np.zeros((N,M))
     for j in range(0,M):
        D2[j,j] = d[j]
     
     '''' create matrices needed to manufacture the low rank matrix'''
     A = np.random.rand(N,N)
     Q1, R = la.qr(A)
     test = np.matmul(Q1,R)
     A =    np.random.rand(M,M)
     Q2,R = la.qr(A)
     
     test = np.matmul(Q2,R)
     
     B = np.matmul(Q1,D2)
     B = np.matmul(B,Q2)
     return B     
     

if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()       
