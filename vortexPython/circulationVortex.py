# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 18:03:34 2015

@author: closadalastra
"""

import numpy as N
import math as M
from numpy import linalg as LA

def circulationPoint(vel_seed_time):
    """
    function that integrates the velocity of a point with time in order to obtain its circulation.
    ---Imputs---
    vel_seed_time - velocity at each point allong the path of the seed with time
    ---Output---
    circulation - approximate value of the circulation at the point
    ---
    C.Losada de la Lastra 2015
    """

    sx,sy,sz = vel_seed_time
    circulation = M.fsum([LA.norm([sx[i],sy[i],sz[i]]) for i in range(len(sx))])
    
    return circulation
    

#def circulationField(L2,p,circ_p,ijk_p):
#    """
#    function that generates the circulation field and adds the fis
#    ---Imputs---
#    vel_seed_time - velocity at each point allong the path of the seed with time
#    ---Output---
#    circulation - approximate value of the circulation at the point
#    ---
#    C.Losada de la Lastra 2015
#    """

def circulationFieldUpdate(pnts,pnts_ijk,pnts_circ,pnts_vect_swirl,CircF,PosF,L2):
    """
    function that updates the circulation field at a certain point.
    ---Imputs---
    CircF - circulation field which will be updated
    p - point/s where the the circulation needs to be upated. 
        [[x0,y0,z0],...,[xn,yn,zn]]
    ijk_p - indices of the point/s where the circulation is calculated in the field 'Circ'
            [[i0,j0,k0],...,[in,jn,kn]]
    circ_p - approximated value of the circulation at the point/s [circ_0,...,circ_n]
    PosF - position field [X,Y,Z]
    L2 - Lambda2 scalar field
    ---Output---
    CircF_updated - updated version of the circulation field
    ---
    C.Losada de la Lastra 2015
    """

    X,Y,Z = PosF    

    decay = 2 #-> DECAY DISTANCE FOR CIRCULATION <-#
    
    for i in range(len(pnts)):
        circ_p = pnts[i]
        circ_p_ijk = ic,jc,kc = pnts_ijk[i]
        circ_p_swirl_vect = pnts_vect_swirl[i]
        circ_p_mag = pnts_circ[i]
        

#        if CircF[kc][jc][ic] == 0: #initially unassigned circulation 
        CircF[kc][jc][ic] = circ_p_mag
#            #circFupdate_decay(p,ijk_p,circ_p,vect_swilr,CircF,PosF,L2,decay,average=None):
#            CircF = circFupdate_decay(circ_p, circ_p_ijk, circ_p_mag,\
#                                        circ_p_swirl_vect, CircF,PosF,L2,decay)
#        elif M.copysign(circ_p_mag,CircF[kc][jc][ic]) < 0: #opposite sign circulations
#            print 'Warning: Trying to assign opposite sign circulations\
#                    to point...Check '+str([X[kc][jc][ic],Y[kc][jc][ic],\
#                    Z[kc][jc][ic]])+' in position field'
#        else: #equal sign but allready assigned circulation
#            CircF[kc][jc][ic] = N.mean([circ_p_mag,CircF[kc][jc][ic]])
#            print 'Warning: Trying to initially assign more than one value of\
#                    circulation to point...Check hasPointCirculation in \
#                    geoVortex for point'+str([X[kc][jc][ic],Y[kc][jc][ic],\
#                    Z[kc][jc][ic]])+' in position field'
        CircF = circFupdate_decay(circ_p, circ_p_ijk, circ_p_mag,\
                                        circ_p_swirl_vect, CircF,PosF,L2,decay,average=1)

    
    return CircF
    
#    ip,jp,kp = ijk_p
#    CircF[kp][jp][ip] = circ_p_0  
#        
#    ###decide weather or not circulation will have a linear distribution with distance
#    for k in range(len(CircF)):
#        for j in range(len(CircF[0])):
#            for i in range(len(CircF[0][0])):
#                if L2[kp][jp][ip] < 0: #point MUST have negative L2 scalar
#                    if circ_p_0 == 0: #no asigned circulation
#                        CircF[kp][jp][ip] = circ_p #WOULD BE A PROBLEM FOR TWO CLOSE VORTICES
#                    elif M.copysign(circ_p,circ_p_0) < 0: #opposite sign circulations
#                        pass
#                    else: #equal sign but allready assigned circulation
#                        CircF[kp][jp][ip] = N.mean(circ_p,circ_p_0)
                        
                                            
                    
        
#    CircF_updated = CircF
#    return CircF_updated    
 

def circFupdate_decay(p,ijk_p,circ_p,vect_swilr,CircF,PosF,L2,decay,average=None):

    X,Y,Z = PosF
    ip,jp,kp = ijk_p
    #define distances away from point in order to reduce looping allong field
    # x,y,z = +/- decay, +/- decay, +/- decay
    dxpos = ip+decay
    if dxpos > len(CircF[0][0]): #check out of range for max value of each
        dxpos = len(CircF[0][0])
    dxneg = ip-decay
    if dxneg < 0:
        dxneg = 0
    dypos = jp+decay
    if dypos > len(CircF[0]):
        dypos = len(CircF[0])
    dyneg = jp-decay
    if dyneg < 0:
        dyneg = 0
    dzpos = kp+decay
    if dzpos > len(CircF):
        dzpos = len(CircF)
    dzneg = kp-decay
    if dzneg < 0:
        dzneg = 0

    if average == None:
        for k in range(dzneg,dzpos+1):
            for j in range(dyneg,dypos+1):
                for i in range(dxneg,dxpos+1):
                    p_temp = N.array([X[k][j][i],Y[k][j][i],Z[k][j][i]])
                    p_temp_l2 = L2[k][j][i]
                    if p_temp_l2 < 0: #point is inside the core
                        p_temp2p = p_temp-p
                        if (N.dot(vect_swilr,p_temp2p) == 0).all(): #point in swirl plane (carefull for very close & opposite L2 surfaces)
                            CircF[k][j][i] = circ_p
                        else: #point is not in swirl plane
    #                        pass
                            vect_swirl_unit = vect_swilr/LA.norm(vect_swilr)
                            dist2swirlPlane = N.dot(vect_swirl_unit,p_temp2p)
                            dist = 2
                            CircF[k][j][i] = circ_p*(1-(dist2swirlPlane/dist))
                            #circulation decays with distance in unintegrated parts of -ve_L2 
                    else:
                        pass
    else:
        dist = 2
        for k in range(dzneg,dzpos+1):
            for j in range(dyneg,dypos+1):
                for i in range(dxneg,dxpos+1):
                    p_temp = N.array([X[k][j][i],Y[k][j][i],Z[k][j][i]])
                    p_temp_l2 = L2[k][j][i]
                    if p_temp_l2 < 0: #point is inside the core
                        p_temp2p = p_temp-p
                        if (N.dot(vect_swilr,p_temp2p) == 0).all(): #point in swirl plane (carefull for very close & opposite L2 surfaces)
                            if CircF[k][j][i] != 0:                                
                                CircF[k][j][i] = N.mean([circ_p,CircF[j][k][i]])
                            else:
                                CircF[k][j][i] = circ_p
                        else: #point is not in swirl plane
    #                        pass
                            if round(CircF[k][j][i],3) != 0:                                
                                vect_swirl_unit = vect_swilr/LA.norm(vect_swilr)
                                dist2swirlPlane = N.dot(vect_swirl_unit,p_temp2p)
                                CircF[k][j][i] = circ_p*(1-(dist2swirlPlane/dist))
                            else:
                                vect_swirl_unit = vect_swilr/LA.norm(vect_swilr)
                                dist2swirlPlane = N.dot(vect_swirl_unit,p_temp2p)
                                CircF[k][j][i] = N.mean([circ_p*(1-(dist2swirlPlane/dist)),CircF[j][k][i]])
                                #circulation decays with distance in unintegrated parts of -ve_L2 
                    else:
                        pass

    return CircF
    
   
def hasPointCirculation(CircF,ijk_p):
    """
    function that checks weather a point has been assigned a value for circulation or not
    ---Imputs---
    CircF - circulation field
    ijk_p - indices of the point where the circulation will be investigated
    ---Output---
    point_circ - True/False weather or not the point has or not a value of circulation
    ---
    C.Losada de la Lastra 2015
    """
    ip,jp,kp = ijk_p
    point_circ = False
    if CircF[kp][jp][ip] != 0:
        point_circ = True
    else:
        pass
    return point_circ
