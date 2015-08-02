# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 18:59:48 2015

@author: closadalastra

"""
import numpy as N
from mayavi import mlab
import matplotlib.pyplot as plt

def vortexTube_2d():
    """
    function that generates a 2d vortex field and writes it to a vtk file
    provide check for the functions Lambda2 and vortexCore.
    """
    from vortexCore import vortexCore2
    from lambda2 import lambda2_2d
    
    gamma = 10 #strength of the vortex.
    Rcore = 10 #radious of the vortex core.
    s = 15 #size parameter of the field x*y points.
    centers = N.array([[0,0,0]]) #points where the vortex will be located
    
    [Y,X] = N.mgrid[-s:s+1, -s:s+1]
    [V,U] = N.zeros_like([Y,X],dtype=float)
    
    #vortexCore2(X,Y,C,Rcore,gamma):
    [U,V] = vortexCore2(X,Y,centers,Rcore,gamma)
            
    L2 = lambda2_2d(X,Y,U,V)
     
    l2 = plt.imshow(L2,cmap='OrRd',extent=[-s,s+1,-s,s+1])
    l2bar = plt.colorbar() 
    l2bar.set_label('Lambda2') 
    vF = plt.quiver(X,Y,U,V,pivot='middle',color= None)    
    
    return


    

def vortexCylinder(): #VERTICAL STACK OF 2D FIELDS....3D FIELD
    """
    """
    from vortexCore import vortexCore2
    from lambda2 import lambda2_3d
    
    gamma = 1 #strength of the vortex.
    Rcore = 18 #radious of the vortex core.
    s = 20 #size parameter of the field x*y points.
    h = 3 #min is 1
    q = N.array([[0,0,0]]) #point where the vortex will be located
    
    [z,y,x] = N.mgrid[0:1,-s:s+1, -s:s+1]
    [Z,Y,X] = N.mgrid[-h:h+1,-s:s+1, -s:s+1]
    [W,V,U] = N.zeros_like([Z,Y,X],dtype=float)
    for i in range(len(Z)):
        #vortexCore2(X,Y,C,Rcore,gamma):
        [U_temp,V_temp] = vortexCore2(x[0],y[0],q,Rcore,gamma)
        
        U[i]=U_temp
        V[i]=V_temp
    
    L2 = lambda2_3d(X,Y,Z,U,V,W)
    #print L2[3]
    #vF=mlab.quiver3d(Z,Y,X,W,V,U)
    Lamda2 = mlab.contour3d(Z,Y,X,L2,transparent=True)
    mlab.show()
    return

#vortexCylinder()
#vortexTube_2d()     
