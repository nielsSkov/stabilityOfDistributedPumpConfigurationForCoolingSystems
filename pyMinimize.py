#minimize cost function, two methods [python in matlab]

from sympy import *
import numpy as np
from scipy import optimize

#Nelder-Mead Simplex algorithm
def BFGS( q, Q, w, R, r, a, b, R_c):
    """Solve using the Broyden–Fletcher–Goldfarb–Shanno algorithm"""
    
    q_init = [ q[0], q[1], q[2], q[3] ]
    
    def f(q_init):
        
        q = np.matrix( [[ q_init[0] ],
                        [ q_init[1] ],
                        [ q_init[2] ],
                        [ q_init[3] ]] )
        
        R1 = R[0];  R2 = R[1];  R3 = R[2];  R4 = R[3];
        r1 = r[0];  r2 = r[1];  r3 = r[2];  r4 = r[3];
        w1 = w[0];  w2 = w[1];  w3 = w[2];  w4 = w[3];
        a1 = a[0];  a2 = a[1];  a3 = a[2];  a4 = a[3];
        b1 = b[0];  b2 = b[1];  b3 = b[2];  b4 = b[3];
        
        S1 = np.matrix( [[ R_c/b1 + (a1 + r1)/b1 + 2*R1/b1, R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ],
                         [ R_c/b1 + 2*R1/b1,                R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ],
                         [ R_c/b1 + 2*R1/b1,                R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ],
                         [ R_c/b1 + 2*R1/b1,                R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ]] )


        S2 = np.matrix( [[ R_c/b2 + 2*R1/b2,  R_c/b2 + 2*R1/b2,                          R_c/b2 + 2*R1/b2,            R_c/b2 + 2*R1/b2           ],
                         [ R_c/b2 + 2*R1/b2,  R_c/b2 + (a2 + r2)/b2 + 2*R1/b2 + 2*R2/b2, R_c/b2 + 2*R1/b2 + 2*R2/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2 ],
                         [ R_c/b2 + 2*R1/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2,                R_c/b2 + 2*R1/b2 + 2*R2/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2 ],
                         [ R_c/b2 + 2*R1/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2,                R_c/b2 + 2*R1/b2 + 2*R2/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2 ]] )


        S3 = np.matrix( [[ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3,            R_c/b3 + 2*R1/b3,                                     R_c/b3 + 2*R1/b3                     ],
                         [ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,                           R_c/b3 + 2*R1/b3 + 2*R2/b3           ],
                         [ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,  R_c/b3 + (a3 + r3)/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3 ],
                         [ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3,                 R_c/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3 ]] )


        S4 = np.matrix( [[ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4,            R_c/b4 + 2*R1/b4,                      R_c/b4 + 2*R1/b4                                              ],
                         [ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,            R_c/b4 + 2*R1/b4 + 2*R2/b4                                    ],
                         [ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4                          ],
                         [ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4,  R_c/b4 + (a4 + r4)/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4 + 2*R4/b4 ]] )
        
        J_cost = ( q.T*S1*q - w1**2 )**2 + \
                 ( q.T*S2*q - w2**2 )**2 + \
                 ( q.T*S3*q - w3**2 )**2 + \
                 ( q.T*S4*q - w4**2 )**2
        
        #print(J_cost.item())
        
        return J_cost.item()
    
    
    #optResult = optimize.minimize(f, q_init, method="BFGS", options={'gtol': 2e-2, 'disp': False})
    
    qBounds = optimize.Bounds(0, np.inf)

    #Broyden–Fletcher–Goldfarb–Shanno algorithm
    optResult = optimize.minimize( f, q_init, method="L-BFGS-B", bounds=qBounds,            \
                                   options={ 'disp':                 None,                  \
                                             'maxcor':               10,                    \
                                             'ftol':                 2.220446049250313e-09, \
                                             'gtol':                 1e-05,                 \
                                             'eps':                  1e-08,                 \
                                             'maxfun':               15000,                 \
                                             'maxiter':              15000,                 \
                                             'iprint':               -1,                    \
                                             'maxls':                20,                    \
                                             'finite_diff_rel_step': None                   })
    
    q_out     = optResult.x
    error     = f(q_out)
    
    return [ q_out[0], q_out[1], q_out[2], q_out[3] ]


#Nelder-Mead Simplex algorithm
def nelderMead( q, Q, w, R, r, a, b, R_c ):
    """Solve using the Nelder-Mead Simplex algorithm"""
    
    q_init = [ q[0], q[1], q[2], q[3] ]
    
    def f(q_init):
        
        q = np.matrix( [[ q_init[0] ],
                        [ q_init[1] ],
                        [ q_init[2] ],
                        [ q_init[3] ]] )
        
        R1 = R[0];  R2 = R[1];  R3 = R[2];  R4 = R[3];
        r1 = r[0];  r2 = r[1];  r3 = r[2];  r4 = r[3];
        w1 = w[0];  w2 = w[1];  w3 = w[2];  w4 = w[3];
        a1 = a[0];  a2 = a[1];  a3 = a[2];  a4 = a[3];
        b1 = b[0];  b2 = b[1];  b3 = b[2];  b4 = b[3];
        
        S1 = np.matrix( [[ R_c/b1 + (a1 + r1)/b1 + 2*R1/b1, R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ],
                         [ R_c/b1 + 2*R1/b1,                R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ],
                         [ R_c/b1 + 2*R1/b1,                R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ],
                         [ R_c/b1 + 2*R1/b1,                R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1,  R_c/b1 + 2*R1/b1 ]] )


        S2 = np.matrix( [[ R_c/b2 + 2*R1/b2,  R_c/b2 + 2*R1/b2,                          R_c/b2 + 2*R1/b2,            R_c/b2 + 2*R1/b2           ],
                         [ R_c/b2 + 2*R1/b2,  R_c/b2 + (a2 + r2)/b2 + 2*R1/b2 + 2*R2/b2, R_c/b2 + 2*R1/b2 + 2*R2/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2 ],
                         [ R_c/b2 + 2*R1/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2,                R_c/b2 + 2*R1/b2 + 2*R2/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2 ],
                         [ R_c/b2 + 2*R1/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2,                R_c/b2 + 2*R1/b2 + 2*R2/b2,  R_c/b2 + 2*R1/b2 + 2*R2/b2 ]] )


        S3 = np.matrix( [[ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3,            R_c/b3 + 2*R1/b3,                                     R_c/b3 + 2*R1/b3                     ],
                         [ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,                           R_c/b3 + 2*R1/b3 + 2*R2/b3           ],
                         [ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,  R_c/b3 + (a3 + r3)/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3 ],
                         [ R_c/b3 + 2*R1/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3,  R_c/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3,                 R_c/b3 + 2*R1/b3 + 2*R2/b3 + 2*R3/b3 ]] )


        S4 = np.matrix( [[ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4,            R_c/b4 + 2*R1/b4,                      R_c/b4 + 2*R1/b4                                              ],
                         [ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,            R_c/b4 + 2*R1/b4 + 2*R2/b4                                    ],
                         [ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4                          ],
                         [ R_c/b4 + 2*R1/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4,  R_c/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4,  R_c/b4 + (a4 + r4)/b4 + 2*R1/b4 + 2*R2/b4 + 2*R3/b4 + 2*R4/b4 ]] )
        
        J_cost = ( q.T*S1*q - w1**2 )**2 + \
                 ( q.T*S2*q - w2**2 )**2 + \
                 ( q.T*S3*q - w3**2 )**2 + \
                 ( q.T*S4*q - w4**2 )**2
        
        #print(J_cost.item())
        
        return J_cost.item()
    
    optResult = optimize.minimize(f, q_init, method='nelder-mead', options={'xatol': 1e-8, 'disp': False})
    
    q_out     = optResult.x
    error     = f(q_out)
    
    return [ q_out[0], q_out[1], q_out[2], q_out[3] ]


