#!/usr/bin/env python3

from sympy import *
import numpy as np
from scipy import signal
from scipy import linalg

def lqr(A,B,Q,R): #contineous time LQR

    #solve ricatti equation
    X = np.matrix( linalg.solve_continuous_are( A, B, Q, R ) )
     
    #compute LQR gain
    K = np.matrix( linalg.inv(R)*(B.T*X) )
     
    eigVals, eigVecs = linalg.eig( A - B*K )
     
    return K, X, eigVals
 
def dlqr(A,B,Q,R): #discrete time LQR
 
    #solve ricatti equation
    X = np.matrix( linalg.solve_discrete_are( A, B, Q, R ) )
     
    #compute LQR gain
    K = np.matrix( linalg.inv( B.T*X*B + R )*( B.T*X*A ) )
     
    eigVals, eigVecs = linalg.eig( A - B*K )
     
    return K, X, eigVals


#### SYSTEM ###########################################################################

A = np.array([[ -0.0499546966880376, 0.0249773146037843, 0 ],
              [  2.99020707183973,  -3.98694410887676,   0 ],
              [  0,                  1,                  0 ]])

B = np.array([[ -26.0416314887154 ],
              [  0               ],
              [  0               ]])

P = np.array([ -1, -2, -3 ])


#### POLE PLACEMENT DESIGN ############################################################

#compute gain matrix, K, in 3 ways
K1 = signal.place_poles(A, B, P, method='KNV0')
K2 = signal.place_poles(A, B, P, method='YT')
K3 = signal.place_poles(A, B, P, method='YT', rtol=-1, maxiter=100)

#print the found gain matrices to terminal
print(K1.gain_matrix)
print(K2.gain_matrix)
print(K3.gain_matrix)


#### LINEAR QUADRATIC REGULATOR (LQR) DESIGN ##########################################


#desired max for system states
x1Max = 2;  # water temperature difference from equilibrium, theta-theta*
x2Max = 1;  # temperature difference from equilibrium, T-T*
x3Max = 30;  # temperature integral state, integral( T-T* )

#desired max actuation (max water flow, q) 
uMax = 100;

Q = np.array([[ 1/(x1Max**2), 0,            0            ],
              [ 0,            1/(x2Max**2), 0            ],
              [ 0,            0,            1/(x3Max**2) ]])

R = np.array([[ uMax ]])

K_lqr, X_lqr, eigVals_lqr = lqr(A,B,Q,R)

#print the found gain martix
print("LQR Gain Matrix")
print( K_lqr )
print("LQR Eigenvalues")  # The system is stable if and only if
print( eigVals_lqr[0] )   # all of the eigenvalues have
print( eigVals_lqr[1] )   # strictly negative real parts
print( eigVals_lqr[2] )













