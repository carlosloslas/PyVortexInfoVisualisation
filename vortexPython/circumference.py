# -*- coding: utf-8 -*-
"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
""" 
import numpy as N
import math as M
#from mayavi import mlab
from numpy import linalg as LA
    
def circumference(c,r,Nvec,np):
    """
    Function that generates the coordinates of a circunference as a 
    3xnumberOfPoints matrix.
    ---Inputs---
    c - center of circumfarence.
    r - radius of circumference.
    Nvec - vector normal to the circumference, used to rotate the circunfarence (vortex ring)
    np - number of points defining the circumference
    ---
    C.Losada de la Lastra 2014
    
    NOTE: For the purpose of this project the matrix 'C' contains np+1 columns.
    This is to take into account the need to close the circunference when 
    plotting the circunference or calculating the induced velocity due to the vortex ring.
    """
    from rank_nullspace import nullspace #.pyc file where nullspace is defined

    #range of angles corresponding to the points on the circumference
    theta = N.array([i*(2*M.pi/np) for i in range(np)])
        
    thet_cos = N.array([M.cos(theta[i]) for i in range(len(theta))])
    theta_cos = N.array([thet_cos,])
    thet_sin = N.array([M.sin(theta[i]) for i in range(len(theta))])
    theta_sin = N.array([thet_sin,])
    #rotation of the circunference with respect to the vector normal to the ring    
    rot = N.array(nullspace(Nvec))#rotation vecors y & z  
    rot_y = N.array([rot[:,0],]).transpose()
    rot_z = N.array([rot[:,1],]).transpose()
    
    #rotation of the x,y,z coordinates of the cirumference
    C_origin =  N.array([c,]*(np)).transpose()
    C_y = r*N.dot(rot_y,theta_cos)
    C_z = r*N.dot(rot_z,theta_sin) 

    C = C_origin + C_y + C_z
    #points in the middle of the segments forming the circumference
    #where the velocity will be induced from (more realistic)    
    C_iVel = N.zeros_like(C)
    for i in range(len(C_y[0])):
        pnt = N.array([C[0][i],C[1][i],C[2][i]])
        pnt0 = N.array([C[0][i-1],C[1][i-1],C[2][i-1]])
        vect = pnt0-pnt
        vect_u = vect/LA.norm(vect)
        
        ivel_pnt = pnt + vect_u*(LA.norm(vect)/2)
        
        C_iVel[0][i]=ivel_pnt[0]        
        C_iVel[1][i]=ivel_pnt[1]        
        C_iVel[2][i]=ivel_pnt[2]        
   
    return N.array([C, C_iVel])
   