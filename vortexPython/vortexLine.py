# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 21:53:00 2015

@author: closadalastra
"""
def vortexLine3(PosF,C,TVR,VR_iVel,gamma):
    """
    function that generates the induced velocity of a vortex line
    the vortex line has an infinitely small core radious
    (developed for vortex rings with no core width)
    ---Inputs---
    PosF - Position field 3d arrays. [X,Y,Z]
    C - Array with the centers of the rings. [[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]]
    TVR - Array with the tangent vectors allong the vortex ring. 
    VR_iVel - Array of the midpoints allong the vortex ring segments, where the induced velocity is calculated
    gamma - strength of the vortex ring.
    ---Outputs---
    [U,V,W] - Induced veloctity field due to the vortex rings.
    ---
    C.Losada de la Lastra 2015
    """
    import numpy as N
    from biotSavart import biotSavart 
    
    X,Y,Z = PosF    
    [W,V,U] = N.zeros_like([Z,Y,X],dtype=float)
    
    #loop thought points with i j k coordinates
    for k in range(len(Z)):
        print(str(k)+' out of '+str(len(Z)))
        for j in range(len(Z[0])-1):
            for i in range(len(Z[0][0])):
                p_temp = N.array([X[k][j][i],Y[k][j][i],Z[k][j][i]])
                #loop though number of rings
                for n in range(C[0][0]):
                    u_temp = 0.
                    v_temp = 0.
                    w_temp = 0.
                    #loop through each ring and obtain the induced velocity.
                    for m in range(len(TVR[0][0])):
                        #tangent vector
                        tvr_temp = N.array([TVR[n][0][m], TVR[n][1][m], TVR[n][2][m]])
                        #point where the velocity is calculated about
                        vr_temp = N.array([VR_iVel[n][0][m],VR_iVel[n][1][m],VR_iVel[n][2][m]])
                        #---biotSavart(tv,q,p, gamma):---#
                        vInd_temp = biotSavart(tvr_temp,vr_temp,p_temp,gamma)
                        
                        u_temp = u_temp + vInd_temp[0]
                        v_temp = v_temp + vInd_temp[1]
                        w_temp = w_temp + vInd_temp[2]
                                      
                    U[k][j][i] = u_temp +U[k][j][i]
                    V[k][j][i] = v_temp +V[k][j][i]
                    W[k][j][i] = v_temp +W[k][j][i]
    return [U,V,W]
                    