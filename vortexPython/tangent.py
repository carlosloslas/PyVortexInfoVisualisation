# -*- coding: utf-8 -*-
"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
"""
import numpy as N

def tangent(C):
    """
    Function that generates the tangential vectors allong a closed curve 'C'
    defined by n number of points.(the start and the end point should be the
    same but they should both be cointained in the 3xn matrix 'C')
    ---Inputs---
    C - Matrix with the points of the closed curve.
    ---Method---
    Obtain tangent vectors:
    P[a,b,c] Q[d,e,f] => PQ[d-a,e-b,c-f]
    ---Output---
    T - tangent vectors allong the curve 'C'
    ---
    C.Losada de la Lastra 2014
       
    NOTE: For the purpose of this project the matrix 'T' contains the initial 
    vector twice -at the begging and at the end-. This is to take into account 
    the need to close the vortex line when calculating the induced velocity.
    """
    t_x=N.array([C[0,i]-C[0,i-1] for i in range(len(C[0]))])
    t_y=N.array([C[1,i]-C[1,i-1] for i in range(len(C[1]))])
    t_z=N.array([C[2,i]-C[2,i-1] for i in range(len(C[2]))])
    
    T = N.vstack((t_x,t_y,t_z))
    
    return T 