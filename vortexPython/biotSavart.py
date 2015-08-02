"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
"""
import numpy as N
from numpy import linalg as LA
import math as M    

def biotSavart(tv,q,p, gamma):
    """
    Function that computes applying the Biot-Savart Law  the velocity induced
    by a vortex line.
    ---Inputs---
    tv - tangent vector to the vortex line
    q - point allong the vortex line where 'tv' is located
    p - point in the field where the induced velocity will be calculated.
    gamma - strength of the vortex line.
    ---
    C.Losada de la Lastra 2014
    """  
    qp = p-q 
    dist = LA.norm(qp)
    qp_unit = pq/dist
    if dist==0.:
        print("Warning:Points p and q are the same point")#debugging & information purposes
        print 'p ='+str(p)
        print 'q ='+str(q)
        return [0,0,0] #Since there will be no induced velocity.
    
    tv_unit = tv/LA.norm(tv)
    #calculating & normalising velocity induced vector.
    vInd_vect = N.cross(qp_unit,tv_unit)
    norm_vInd_vect = LA.norm(vInd_vect)
    if norm_vInd_vect==0.: #catching nan error. ie [0,0,0]/0. = [nan,nan,nan]
        print("Warning:vel induced vector has a length of 0 between p and q")#debugging & information purposes
        print 'p ='+str(p)
        print 'q ='+str(q)
        return [0,0,0] #Since there will be no induced velocity.
    vInd_unit = vInd_vect/norm_vInd_vect
    
    #Applying Biot-Savart Law
    vInd= (gamma/(4.*M.pi*(dist)**2))*vInd_unit #vInd is unit vector so **3 becomes **2.
    
    return vInd