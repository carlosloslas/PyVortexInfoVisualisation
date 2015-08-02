# -*- coding: utf-8 -*-
"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
"""
import numpy as N
from numpy import linalg as LA
import math as M    

############################
#-VORTEX LINE WITH CORE 3D-#
############################
def vortexCore3(PosF,C,TVR,VR_iVel,Rcore,gamma):
    """
    Function that generates the induced 3D velocity field due to a vortex ring
    ---Inputs---
    PosF - Position field 3d arrays. [X,Y,Z]
    C - Array with the centers of the rings. [[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]]
    TVR - array of the tangent vectors allong the rings
    VR_iVel - array of the mid-segment point allong the vortex rings
    Rcore - radius of the core of the vortex rings
    gamma - strength of the vortex rings
    ---Outputs---
    [U,V,W] - Induced veloctity field due to the vortex rings.
    ---
    C.Losada de la Lastra 2015
    """
    from biotSavart import biotSavart
    
    [X,Y,Z] = PosF
    [W,V,U] = N.zeros_like([Z,Y,X],dtype=float)    
    #loop though number of rings
    for n in range(C[0][0]):
        print(str(n)+' out of '+str(C[0][0]-1))
        #loop thought points with i j k coordinates
        for k in range(len(Z)):
            
            for j in range(len(Z[0])):
                for i in range(len(Z[0][0])):
                    p_temp = N.array([X[k][j][i],Y[k][j][i],Z[k][j][i]])
                    
                    u_temp = 0.
                    v_temp = 0.
                    w_temp = 0.
                
                    #loop through each ring and obtain the induced velocity.
                    for m in range(len(TVR[0][0])):
                        #tangent vector
                        tvr_temp = N.array([TVR[n][0][m], TVR[n][1][m], TVR[n][2][m]])
                        #point where the velocity is calculated about
                        vr_temp = N.array([VR_iVel[n][0][m],VR_iVel[n][1][m],VR_iVel[n][2][m]])
                        dist = LA.norm(vr_temp - p_temp)
                        if dist > Rcore:                            
                            #---biotSavart(tv,q,p, gamma):---#
                            #[u_temp,v_temp,w_temp] = biotSavart(tvr_temp,vr_temp,p_temp,gamma)
                            
                            #---biot-savart coupled to rankine function---#
                            #biotSavartRankine(tv,q,p, gamma,Rcore):
                            [u_temp,v_temp,w_temp] = biotSavartRankine(tvr_temp,vr_temp,p_temp, gamma,Rcore)
                        
                        if dist <= Rcore:
                            [u_temp,v_temp,w_temp] = rankineVortex(tvr_temp,vr_temp,p_temp,gamma,Rcore)
                        
                        U[k][j][i] = U[k][j][i] + u_temp
                        V[k][j][i] = V[k][j][i] + v_temp
                        W[k][j][i] = W[k][j][i] + w_temp
    return [U,V,W]

##############################
#  VORTEX LINE WITH CORE 2D  #
##############################
def vortexCore2(X,Y,C,Rcore,gamma):
    """
    Function that generates the induced 3D velocity field due to a vortex ring
    ---Inputs---
    X - X-coordinate matrix
    Y - Y-coordinate matrix
    C - array of the centers of the rings
    Rcore - rafious of the core of the vortex rings
    gamma - strength of the vortex rings
    ---Outputs---
    [U,V] - matrices of the 2d induced velocity field
    ---
    C.Losada de la Lastra 2015
    """
#    from biotSavart import biotSavart
    
    #since it is a 2d field on the xy plane the vector ...
    # on the direction of the vortex line is:
    n_vect = [0,0,1]
    [V,U] = N.zeros_like([Y,X],dtype=float)
    
    #loop though number of vortex centers
    for n in range(len(C)):
        #loop thought points with j k coordinates
        for j in range(len(Y)):
            for i in range(len(Y[0])):
                p_temp = N.array([X[j][i],Y[j][i],0])
                u_temp = 0.
                v_temp = 0.
                dist=LA.norm(p_temp-C[n])
                
                if dist <= Rcore:
                    #rankineVortex(tv,q,p,gamma,Rcore):
                    [u_temp,v_temp,w_temp] = rankineVortex(n_vect,C[n],p_temp,gamma,Rcore)
                else:
                    #biot-savart coupled to the rankine vortex.
                    #biotSavartRankine(tv,q,p, gamma,Rcore):
                    [u_temp,v_temp,w_temp] = biotSavartRankine(n_vect,C[n],p_temp,gamma,Rcore)
                    
                U[j][i] = U[j][i] + u_temp
                V[j][i] = V[j][i] + v_temp
    
    return [U,V]

###############################
#  RANKINE VORTEX DEFINITION  #
###############################
def rankineVortex(tv,q,p,gamma,Rcore):
    """
    function that calculates the induced velocity field inside the core of a rankine vortex
    ---Inputs---
    tv - tangent vector allong vortex line.
    q - point at allong the vortex line @ the tangent vector.
    p - point where the induced velocity is analysed. 
    gamma - strength of the vortex line.
    Rcore - radious of the core of the vortex line.
    ---Outputs---
    vInd - [u,v,w]    
    ---
    C.Losada de la Lastra 2015
    """
    pq = q-p 
    dist = LA.norm(pq)
    pq_unit = pq/dist
    if dist==0.:
        print("Warning:Points p and q are the same point")#debugging & information purposes
        print 'p ='+str(p)
        print 'q ='+str(q)
        return [0,0,0] #Since there will be no induced velocity.
    tv_unit = tv/LA.norm(tv)
    
    #calculating & normalising velocity induced vector.
    vInd_vect = N.cross(pq_unit,tv_unit)

    norm_vInd_vect = LA.norm(vInd_vect)
    if norm_vInd_vect==0.: #catching nan error. ie [0,0,0]/0. = [nan,nan,nan]
        print("Warning:vel induced vector has a length of 0 between p and q")#debugging & information purposes
        print 'p ='+str(p)
        print 'q ='+str(q)
        return [0,0,0] #Since there will be no induced velocity.
    vInd_unit = vInd_vect/norm_vInd_vect
    
    #magnitude * unit vector in direction of induced velocity.
    vInd = (gamma/(4.*M.pi)*(dist/Rcore))*vInd_unit 

    return vInd
    
###########################################
#  BIOT-SAVART COUPLED TO RANKINE VORTEX  #
###########################################
def biotSavartRankine(tv,q,p, gamma,Rcore):
    """
    function that couples the biot-savart function to the rankine function
    ---Inputs---
    tv - tangent vector allong vortex line.
    q - point at allong the vortex line @ the tangent vector.
    p - point where the induced velocity is analysed. 
    gamma - strength of the vortex line.
    Rcore - radious of the core of the vortex line.
    ---Outputs---
    vInd - [u,v,w]    
    ---
    C.Losada de la Lastra 2015
    """       
    pq = q-p 

    dist = LA.norm(pq) #-(Rcore*0.9) 
    #taken away the radious of the core in order to couple both flow regimes.
    
    pq_unit = pq/dist
    if dist==0.: #catching nan error. ie [0,0,0]/0. = [nan,nan,nan]
        print("Warning:Points p and q are the same point")#debugging & information purposes
        print 'p ='+str(p)
        print 'q ='+str(q)
        return [0,0,0] #Since there will be no induced velocity.
    
    tv_unit = tv/LA.norm(tv)
    
   #generating and normalysing the induced velocity vector.
    vInd_vect = N.cross(pq_unit,tv_unit)
    norm_vInd_vect = LA.norm(vInd_vect)
    if norm_vInd_vect==0.: #catching nan error. ie [0,0,0]/0. = [nan,nan,nan]
        print("Warning:vel induced vector has a length of 0 between p and q")#debugging & information purposes
        print 'p ='+str(p)
        print 'q ='+str(q)
        return [0,0,0] #Since there will be no induced velocity.        
    vInd_unit = vInd_vect/norm_vInd_vect    

    #Applying Biot-Savart Law
    vInd= (gamma/(4.*M.pi))*(Rcore/((dist)**2))*(vInd_unit) #vInd is unit vector so **3 becomes **2.
    

    return vInd