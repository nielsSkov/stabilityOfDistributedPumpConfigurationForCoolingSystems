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

# A = np.array([[ -0.0499546966880376, 0.0249773146037843, 0 ],
              # [  2.99020707183973,  -3.98694410887676,   0 ],
              # [  0,                  1,                  0 ]])

# B = np.array([[ -26.0416314887154 ],
              # [  0               ],
              # [  0               ]])

P = np.array([ -1, -2, -3 ])

A1 = np.array([[ -0.112592946471008,  0.025020753411791, 0 ],
               [  2.990207071839725, -3.986941639740960, 0 ],
               [                  0,                  1, 0 ]])
               
B1 = np.array([[ -5.797124303542550 ],
               [                  0 ],
               [                  0 ]])


A2 = np.array([[ -0.112733336529734,  0.025051873751855, 0 ],
               [  3.018959062915107, -4.025278507359552, 0 ],
               [                   0,                  1, 0 ]])

B2 = np.array([[ -9.950257172471055 ],
               [                  0 ],
               [                  0 ]])


A3 = np.array([[ -0.112592946471008,  0.025020753411791, 0 ],
               [  2.990207071839725, -3.986941639740960, 0 ],
               [                  0,                  1, 0 ]])

B3 = np.array([[ -11.594248607085101 ],
               [                   0 ],
               [                   0 ]])


A4 = np.array([[ -0.112398447127517,  0.024977314603784, 0 ],
               [  2.990207071839725, -3.986944108876762, 0 ],
               [                  0,                  1, 0 ]])

B4 = np.array([[ -6.944411611690053 ],
               [                  0 ],
               [                  0 ]])


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

A_all = ( A1, A2, A3, A4 );
B_all = ( B1, B2, B3, B4 );

for A, B, i in zip(A_all, B_all, range(1,5)):

    K_lqr, X_lqr, eigVals_lqr = lqr(A,B,Q,R)
    
    print("\n-------------------------------------------------------\n")
    print("Results from system: A" + str(i) + ", B" + str(i))
    
    #print the found gain martix
    print("LQR Gain Matrix")
    print( K_lqr )
    print("\nLQR Eigenvalues")  # The system is stable if and only if
    print( eigVals_lqr[0] )     # all of the eigenvalues have
    print( eigVals_lqr[1] )     # strictly negative real parts
    print( eigVals_lqr[2] )


#### POLE PLACEMENT DESIGN ############################################################

#compute gain matrix, K, in 3 ways
# K1 = signal.place_poles(A1, B1, P, method='KNV0')
# K2 = signal.place_poles(A1, B1, P, method='YT')
# K3 = signal.place_poles(A1, B1, P, method='YT', rtol=-1, maxiter=100)

#print the found gain matrices to terminal
# print(K1.gain_matrix)
# print(K2.gain_matrix)
# print(K3.gain_matrix)
















