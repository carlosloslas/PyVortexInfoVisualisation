# -*- coding: utf-8 -*-
"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
"""
import numpy as N
from numpy import linalg as LA

###################
#  LAMBDA2 IN 3D  #
###################
def lambda2_3d(X,Y,Z,U,V,W):
    """
    Function that returns the Lambda 2 scalar at a point in a 3d field.    
    ---Imputs---
    X - x-coordinates of the field as a 3d matrix.
    Y - y-coordinates of the field as a 3d matrix.
    Z - z-coordinates of the field as a 3d matrix.
    U - u components at the points in the velocity field as a 3d matrix
    V - v components at the points in the velocity field as a 3d matrix
    W - w components at the points in the velocity field as a 3d matrix
    ---
    C.Losada de la Lastra 2015
    """    
    L2 = N.zeros_like(X,dtype=float) 
    #loop through velocity field and calculate the S and Ω matrices.
    #-2 and +1 to avoid 'out of range' error.
    for k in range(len(Z)-2):
        k=k+1
        for j in range(len(Z[0])-2):
            j=j+1
            for i in range(len(Z[0][0])-2):
                i=i+1

                #calculate gradients of velocity field
                ux = [U[k][j][i-1],U[k][j][i+1]], [X[k][j][i-1],X[k][j][i+1]]
                uy = [U[k][j-1][i],U[k][j+1][i]], [Y[k][j-1][i],Y[k][j+1][i]]
                uz = [U[k-1][j][i],U[k+1][j][i]], [Z[k-1][j][i],Z[k+1][j][i]]
                vx = [V[k][j][i-1],V[k][j][i+1]], [X[k][j][i-1],X[k][j][i+1]]
                vy = [V[k][j-1][i],V[k][j+1][i]], [Y[k][j-1][i],Y[k][j+1][i]]
                vz = [V[k-1][j][i],V[k+1][j][i]], [Z[k-1][j][i],Z[k+1][j][i]]
                wx = [W[k][j][i-1],W[k][j][i+1]], [X[k][j][i-1],X[k][j][i+1]]
                wy = [W[k][j-1][i],W[k][j+1][i]], [Y[k][j-1][i],Y[k][j+1][i]]
                wz = [W[k-1][j][i],W[k+1][j][i]], [Z[k-1][j][i],Z[k+1][j][i]]
                                
                #definig the S and Ω matrix elements and matrices
                s11 = 1./2*(dd(ux)+dd(ux))
                s12 = 1./2*(dd(uy)+dd(vx))
                s13 = 1./2*(dd(uz)+dd(wx))
                s21 = 1./2*(dd(vx)+dd(uy))
                s22 = 1./2*(dd(vy)+dd(vy))
                s23 = 1./2*(dd(vz)+dd(wy))
                s31 = 1./2*(dd(wx)+dd(uz))
                s32 = 1./2*(dd(wy)+dd(vz))
                s33 = 1./2*(dd(wz)+dd(wz))
                
                S = N.array([[s11,s12,s13],[s21,s22,s23],[s31,s32,s33]])
                
                om11 = 1./2*(dd(ux)-dd(ux))
                om12 = 1./2*(dd(uy)-dd(vx))
                om13 = 1./2*(dd(uz)-dd(wx))
                om21 = 1./2*(dd(vx)-dd(uy))
                om22 = 1./2*(dd(vy)-dd(vy))
                om23 = 1./2*(dd(vz)-dd(wy))
                om31 = 1./2*(dd(wx)-dd(uz))
                om32 = 1./2*(dd(wy)-dd(vz))
                om33 = 1./2*(dd(wz)-dd(wz))
                
                OM = N.array([[om11,om12,om13],[om21,om22,om23],[om31,om32,om33]])
               
                L = N.dot(S,S)+N.dot(OM,OM)      
                #eigenvalues of l and L2 value
                eig=LA.eig(L)
                Eig=N.sort(eig[0])
                l2=Eig[1]                
                L2[k][j][i] = l2

    return L2 

###################
#  LAMBDA2 IN 2D  #
###################
    
def lambda2_2d(X,Y,U,V):
    """
    Function that returns the Lambda 2 scalar at a point in a 2d field.    
    ---Imputs---
    X - x-coordinates of the field as a 3d matrix.
    Y - y-coordinates of the field as a 3d matrix.
    U - u components at the points in the velocity field as a 3d matrix
    V - v components at the points in the velocity field as a 3d matrix
    ---
    C.Losada de la Lastra
    """    
    L2 = N.zeros_like(X,dtype=float) 
    #loop through velocity field and calculate the S and Ω matrices.
    #-2 and +1 to avoid 'out of range' error.
    for j in range(len(X)-2):
        j=j+1
        for i in range(len(X[0])-2):
            i=i+1
            #calculate gradients of velocity field
            ux = [U[j][i-1],U[j][i+1]], [X[j][i-1],X[j][i+1]]
            uy = [U[j-1][i],U[j+1][i]], [Y[j-1][i],Y[j+1][i]]
            vx = [V[j][i-1],V[j][i+1]], [X[j][i-1],X[j][i+1]]
            vy = [V[j-1][i],V[j+1][i]], [Y[j-1][i],Y[j+1][i]]
           
            #definig the S and Ω matrix elements and matrices
            s11 = 1./2*(dd(ux)+dd(ux))
            s12 = 1./2*(dd(uy)+dd(vx))
            s21 = 1./2*(dd(vx)+dd(uy))
            s22 = 1./2*(dd(vy)+dd(vy))
            
            S = N.array([[s11,s12],[s21,s22]])
            
            om11 = 1./2*(dd(ux)-dd(ux))
            om12 = 1./2*(dd(uy)-dd(vx))
            om21 = 1./2*(dd(vx)-dd(uy))
            om22 = 1./2*(dd(vy)-dd(vy))
            
            OM = N.array([[om11,om12],[om21,om22]])
            
            L = N.dot(S,S)+N.dot(OM,OM)
        
            #eigenvalues of l and L2 value
            eig=LA.eig(L)
            Eig=N.sort(eig[0])
            l2=Eig[1]
            
            L2[j][i] = l2

    return L2 
    
####################################
#  LAMBDA2 SWIRL PLANE ON 3D FIELD #
####################################

def lambda2_3d_swirl(PosF,VelF):
    """
    function that calculates the eigenvectors corresponding to the lambda2 eigenvalues
    it should be run instead of the lambda2_3d function, since it recalculates L2
    
    """
    X,Y,Z = PosF
    U,V,W = VelF

    L2 = N.zeros_like(X,dtype=float)
    L2x = N.zeros_like(X,dtype=float)
    L2y = N.zeros_like(X,dtype=float)
    L2z = N.zeros_like(X,dtype=float)
    L2u = N.zeros_like(X,dtype=float)
    L2v = N.zeros_like(X,dtype=float)
    L2w = N.zeros_like(X,dtype=float)
    
    #loop through velocity field and calculate the S and Ω matrices.
    #-2 and +1 to avoid 'out of range' error.
    for k in range(len(Z)-2):
        k=k+1
        for j in range(len(Z[0])-2):
            j=j+1
            for i in range(len(Z[0][0])-2):
                i=i+1

                #calculate gradients of velocity field
                ux = [U[k][j][i-1],U[k][j][i+1]], [X[k][j][i-1],X[k][j][i+1]]
                uy = [U[k][j-1][i],U[k][j+1][i]], [Y[k][j-1][i],Y[k][j+1][i]]
                uz = [U[k-1][j][i],U[k+1][j][i]], [Z[k-1][j][i],Z[k+1][j][i]]
                vx = [V[k][j][i-1],V[k][j][i+1]], [X[k][j][i-1],X[k][j][i+1]]
                vy = [V[k][j-1][i],V[k][j+1][i]], [Y[k][j-1][i],Y[k][j+1][i]]
                vz = [V[k-1][j][i],V[k+1][j][i]], [Z[k-1][j][i],Z[k+1][j][i]]
                wx = [W[k][j][i-1],W[k][j][i+1]], [X[k][j][i-1],X[k][j][i+1]]
                wy = [W[k][j-1][i],W[k][j+1][i]], [Y[k][j-1][i],Y[k][j+1][i]]
                wz = [W[k-1][j][i],W[k+1][j][i]], [Z[k-1][j][i],Z[k+1][j][i]]
                                
                #definig the S and Ω matrix elements and matrices
                s11 = 1./2*(dd(ux)+dd(ux))
                s12 = 1./2*(dd(uy)+dd(vx))
                s13 = 1./2*(dd(uz)+dd(wx))
                s21 = 1./2*(dd(vx)+dd(uy))
                s22 = 1./2*(dd(vy)+dd(vy))
                s23 = 1./2*(dd(vz)+dd(wy))
                s31 = 1./2*(dd(wx)+dd(uz))
                s32 = 1./2*(dd(wy)+dd(vz))
                s33 = 1./2*(dd(wz)+dd(wz))
                
                S = N.array([[s11,s12,s13],[s21,s22,s23],[s31,s32,s33]])
                
                om11 = 1./2*(dd(ux)-dd(ux))
                om12 = 1./2*(dd(uy)-dd(vx))
                om13 = 1./2*(dd(uz)-dd(wx))
                om21 = 1./2*(dd(vx)-dd(uy))
                om22 = 1./2*(dd(vy)-dd(vy))
                om23 = 1./2*(dd(vz)-dd(wy))
                om31 = 1./2*(dd(wx)-dd(uz))
                om32 = 1./2*(dd(wy)-dd(vz))
                om33 = 1./2*(dd(wz)-dd(wz))
                
                OM = N.array([[om11,om12,om13],[om21,om22,om23],[om31,om32,om33]])
               
                L = N.dot(S,S)+N.dot(OM,OM)      
                #eigenvalues of l and L2 value
                eig = LA.eig(L)
                Eig = N.sort(eig[0])
                l2_temp = Eig[1]                
                L2[k][j][i] = l2_temp     
                                                
                #eigenvector corresponding to l2 eigenvalue when inside core
                if l2_temp < 0:
                    for i in range(len(eig[0])):
                        if eig[0][i] == l2_temp:
                            eigVec_l2 = eig[1][i]
                            eigVec_l2_unit = eigVec_l2/LA.norm(eigVec_l2)
#                            L2x[k][j][i] = X[k][j][i]
#                            L2y[k][j][i] = Y[k][j][i]
#                            L2z[k][j][i] = Z[k][j][i]
                            L2u[k][j][i] = eigVec_l2_unit[0]
                            L2v[k][j][i] = eigVec_l2_unit[1]
                            L2w[k][j][i] = eigVec_l2_unit[2]
                        else:
                            pass
                else:
                    pass
                
                
    
    L2Vec = [L2x,L2y,L2z,L2u,L2v,L2w]

    return [L2,L2Vec]
    
    
    
    

###---util functions---###
#1.velocty gradient: dd(ux):
#2.local maximum L2 values: L2max(PosF,L2):
#3.local minimum L2 values : L2min(PosF,L2):
#4.croping lambda2Min field points: cropL2min(dataL2min):

def dd(ux):
    """
    function that calculates the gradient of 'u' in the direction 'x'
    ---Inputs---    
    ux - array formed by u and x.
    u,x = ux
    u - components of the veocity in the direction interested in.
    N.array([ux-1,ux+1])
    x - coordinates of the points in the direction interested in.
    N.array([px-1,px+1])
    ---
    C.Losada de la Lastra
    """
    u,x = ux
    if x[1]==x[0]:#no change in position in said direction, will cause the gracient to be +/-infinite.
        print 'Velocity graident equals zero for'+str(ux)
        return 0.        
    #final - iniital stage of velocity component over the equivatent in displacement  
    return ((u[1]-u[0])/(x[1]-x[0]))

def L2max(PosF,L2):
    """
    function that locates the points allong the edge of the vortex cores,
    where Lambda2 is a local maximum.
    ---Inputs---
    PosF - matrix of the position field [X,Y,Z]
    L2 - matrix of the lambda2 scalar
    ---Outputs---
    [[L2_1,L2_2,...,L2_n]
    [x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]
    [i1,j1,k1],[i2,j2,k2],...,[in,jn,kn]] - list of the value of lamda2, 
    the location of  points on the edge of the vortex core,
    and their index in the position field.
    ---
    C.Losada de la Lastra 2015
    """
    data = [[] for i in range(3)]
    
    for k in range(len(L2)-2):#avoid running out of index
        k = k+1
        for j in range(len(L2[0])-2):
            j = j+1
            for i in range(len(L2[0][0])-2):
                i = i+1
                testVal = L2[k][j][i]
                if testVal > L2[k-1][j][i] and testVal > L2[k+1][j][i]:
                    if testVal > L2[k][j-1][i] and testVal > L2[k][j+1][i]:
                        if testVal > L2[k][j][i-1] and testVal > L2[k][j][i+1]:
                            #append L2, coordinates, and indices to data list
                            data[0].append(testVal)
                            data[1].append([PosF[0][k][j][i],PosF[1][k][j][i],PosF[2][k][j][i]])
                            data[2].append([i,j,k])
    
    return data
    
def L2min(PosF,L2):
    """
    function that locates the points on the vortex cores,
    where Lambda2 is a minimum.
    ---Inputs---
    PosF - matrix of the position field [X,Y,Z]
    L2 - matrix of the lambda2 scalar
    prox - optional prarmeter that indicates the distance from each 
    L2min point to the next one in order to reduce the data output from the fucntion
    ---Outputs---
    [[L2_1,L2_2,...,L2_n]
    [x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]
    [i1,j1,k1],[i2,j2,k2],...,[in,jn,kn]] - list of the value of lamda2, 
    the location of  points on the edge of the vortex core,
    and their index in the position field.
    ---
    C.Losada de la Lastra 2015
    """
    data = [[] for i in range(3)]
    [X,Y,Z] = PosF
    
    
    for k in range(len(L2)-2):#avoid running out of index
        k = k+1
        for j in range(len(L2[0])-2):
            j = j+1
            for i in range(len(L2[0][0])-2):
                i = i+1
                testVal = L2[k][j][i]
                if testVal < L2[k-1][j][i] and testVal < L2[k+1][j][i]:
                    if testVal < L2[k][j-1][i] and testVal < L2[k][j+1][i]:
                        if testVal < L2[k][j][i-1] and testVal < L2[k][j][i+1]:
                            #append L2, coordinates, and indices to data list
                            data[0].append(testVal)
                            data[1].append([X[k][j][i],Y[k][j][i],Z[k][j][i]])
                            data[2].append([i,j,k])
                                
    return data

def cropL2min(dataL2min,PosF,method='negData'):
    """
    function that crops the L2min data in order to account for false local L2min points
    ---Inputs---
    dataL2min - array out of L2min function. [[L2_val],[coord],[index]]
    method - way of croping the data
        negData: leaves only the negative L2 points
        dist: makes sure the distance between the points isn't too small
        auto: which calls first the 'negData' and then the 'dist' mehtod
    PosF - in order to arrive at a sensible distance to crop the data.[X,Y,Z]
    (not needed else)
    ---Outputs---
    dataL2min_croped - the reduced set of L2min data points
    ---
    C.Losada de la Lastra 2015
    """
    dataL2min_croped = [[] for i in range(3)]     
    
    if method == 'negData':
               
        for i in range(len(dataL2min[0])):
            if round(dataL2min[0][i],4) < 0.0: #if current value of L2 is negative
                #add the three pieces of data to the the croped data list.
                dataL2min_croped[0].append(dataL2min[0][i])
                dataL2min_croped[1].append(dataL2min[1][i])
                dataL2min_croped[2].append(dataL2min[2][i])
            else:
                pass
        return dataL2min_croped
        
    if method == 'dist':
        #come up with a reasonable distance between neighbouring L2min points
    
        crit = 3 * 2 # half of a third of the length of a field, 
        #would need to be revisited for improvements
        
        X,Y,Z = PosF
        Xlen = len(X[0][0])
        Ylen = len(X[0])
        Zlen = len(X)
        xlen = Xlen/crit
        ylen = Ylen/crit
        zlen = Zlen/crit
        xl = round(xlen)
        yl = round(ylen)
        zl = round(zlen)
#        xl=5
#        yl=5
#        zl=5
        
        val,pnt,ijk = dataL2min        
        #loop though data and check if distance with the point before
        #is less than than the one above defined.
        i=1
        #add first point of 'dataL2min'
        dataL2min_croped[0].append(dataL2min[0][i-1])
        dataL2min_croped[1].append(dataL2min[1][i-1])
        dataL2min_croped[2].append(dataL2min[2][i-1])
        #loop through the data and compare the distance between the points 
        #with the distance defined initially in the 'method'
        for n in range(1,len(dataL2min[0])):
            if abs(abs(pnt[i-1][0]) - abs(pnt[i][0])) < xl or \
                    abs(abs(pnt[i-1][1]) - abs(pnt[i][1])) < yl or \
                    abs(abs(pnt[i-1][2]) - abs(pnt[i][2])) < zl:
                
                #check which of the two L2min has smaller L2
                if val[i-1] < val[i]:
                    #point on position 'i-1' has allready been added
                    #remove point from data
                    dataL2min[0].pop(i)
                    dataL2min[1].pop(i)
                    dataL2min[2].pop(i)
                else:
                    #add the point since its smaller
                    dataL2min_croped[0].append(dataL2min[0][i])
                    dataL2min_croped[1].append(dataL2min[1][i])
                    dataL2min_croped[2].append(dataL2min[2][i])
                    #remove point before from data
                    dataL2min[0].pop(i-1)
                    dataL2min[1].pop(i-1)
                    dataL2min[2].pop(i-1)
            else:
                i += 1
                pass
        
        return dataL2min_croped
        
    if method == 'auto':
        dataL2min_croped_0 = cropL2min(dataL2min,PosF,method='negData')
        dataL2min_croped = cropL2min(dataL2min_croped_0,PosF,method='dist')        
                    
        return dataL2min_croped
    
        