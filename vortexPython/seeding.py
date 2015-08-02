# -*- coding: utf-8 -*-
"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
"""
import numpy as N

def seed_time(p0,ijk,PosF,VelF,t):
    """
    function that calculates the evolution in time of a seed placed at a point in space
    ---Inputs---
    p0 - initial location of the seed [x0,y0,z0]
    ijk - indices of the point in the PosF and VelF fields.
    PosF - matrix of the position field [X,Y,Z]
    VelF - matrix of the velocity field [U,V,W]
    t - time in seconds
    ---Outputs---
    [ [[x0,y0,z0],[x1,y1,z1],...,[xn,yn,zn]],
    [[u0,v0,z0],[u1,v1,z1],...,[un,vn,zn]] ] - Matrix of  2 x numberOfTimeSteps, 
    with the position of the seed and the velocity of the seed @ each time step.
    ---
    C.Losada de la Lastra 2015
    """
    #DEFINE TIME STEP
    dt = 0.25 #seconds
    nt = int(t/dt) #number of time steps
    #indices of field matrices
    i0,j0,k0 = ijk
    U,V,W = VelF
    vel_dt0 = U[k0][j0][i0], V[k0][j0][i0], W[k0][j0][i0]      
    
    [PosVec,VelVec] = N.zeros((2,3,nt)) #output data array
    ijk_dt = N.zeros((nt,3)) #index of seed path for every dt
    
    PosVec[0][0], PosVec[1][0], PosVec[2][0] = p0
    VelVec[0][0], VelVec[1][0], VelVec[2][0] = vel_dt0
    ijk_dt[0] = i0, j0, k0
    t_dt = 0
    for i in range(1,nt): #loop for every thime step
        t_dt += dt
        p_temp_dt = N.array([PosVec[0][i-1], PosVec[1][i-1], PosVec[2][i-1]])
        ijk_temp_dt = ijk_dt[i-1]
        
        [PosVec_temp, VelVec_temp] = seed_time_2nAprox(p_temp_dt, ijk_temp_dt, PosF, VelF, t_dt)
#        [PosVec_temp, VelVec_temp] = seed_time_3rdAprox(p_temp_dt, ijk_temp_dt, PosF, VelF, t_dt)
        PosVec[0][i],PosVec[1][i],PosVec[2][i] = PosVec_temp
        VelVec[0][i],VelVec[1][i],VelVec[2][i] = VelVec_temp
        ijk_dt[i] = findPosIjk(PosVec_temp,PosF)
   
    return [PosVec, VelVec]

def seed_time_3rdAprox(p0,ijk0,PosF,VelF,dt):
    """
    function that calculates the evolution in time of a seed 
    placed at a point in space.
    ---Method---
    this fucntion uses a third order aproximation for the new position 
    of the seed.
    ---Inputs---
    p0 - location of the point. [xp,yp,zp]
    ijk0 - indices of the point 'p0' in the fields 'PosF' and 'VelF'
    PosF - matrix of the position field [X,Y,Z]
    VelF - matrix of the velocity field [U,V,W]
    dt - differential of time between timesteps
    ---Outputs---
    [[xi,yi,zi], [ui,vi,wi]] - Array containning the new position and velocity after a dt.
    ---
    C.Losada de la Lastra 2015
    """
    i0,j0,k0 = ijk0
    v_0 = N.array([VelF[0][k0][j0][i0],VelF[1][k0][j0][i0],VelF[2][k0][j0][i0]])
    
    #second order approximaiton of position
    x_0 = N.array(p0)
    x_1 = p0 + dt*v_0
    v_1 = findVel(x_1,PosF,VelF)
    x_2 = x_0 + (dt*0.5*v_0) + (dt*0.5*v_1)
    v_2 = findVel(x_2,PosF,VelF)
    
    #third order approximaiton of position
    x_3 = x_0 + ((1/3.)*dt*v_0) +((1/3.)*dt*v_1) +((1/3.)*dt*v_2)
    v_3 = findVel(x_3,PosF,VelF)
    
    return [x_3,v_3]
    
def seed_time_2nAprox(p0,ijk0,PosF,VelF,dt):
    """
    function that calculates the evolution in time of a seed 
    placed at a point in space.
    ---Method---
    this fucntion uses a second order aproximation for the new position 
    of the seed.
    ---Inputs---
    p0 - location of the point. [xp,yp,zp]
    ijk0 - indices of the point 'p0' in the fields 'PosF' and 'VelF'
    PosF - matrix of the position field [X,Y,Z]
    VelF - matrix of the velocity field [U,V,W]
    dt - differential of time between timesteps
    ---Outputs---
    [[xi,yi,zi], [ui,vi,wi]] - Array containning the new position and velocity after a dt.
    ---
    C.Losada de la Lastra 2015
    """
    i0,j0,k0 = ijk0
    v_0 = N.array([VelF[0][k0][j0][i0],VelF[1][k0][j0][i0],VelF[2][k0][j0][i0]])
    
    #second order approximaiton of position
    x_0 = N.array(p0)
    x_1 = p0 + dt*v_0
    v_1 = findVel(x_1,PosF,VelF)
    x_2 = x_0 + (dt*0.5*v_0) + (dt*0.5*v_1)
    v_2 = findVel(x_2,PosF,VelF)
    
    return [x_2, v_2]
    
def findVel(p,PosF,VelF):
    """
    function that aproximates the velocity at the point p in the field VelF
    ---Inputs---
    p - point where the velocity will be found. [x,y,z]
    PosF - 3d array with the position field. [X,Y,Z]
    VelF - 3d array with the velocity field. [U,V,W]
    ---Outputs---
    Vel - velocity at point p. [u,v,w]
    ---
    C.Losada de la Lastra 2015
    """
    #round coordinates of p
    p_round = [round(p[i]) for i in range(len(p))]
    
    #find indices corresponding to p_rounded
    x_coords = PosF[0][0][0]
    y_coords = [PosF[1][0][i][0] for i in range(len(PosF[1][0]))]
    z_coords = [PosF[2][i][0][0] for i in range(len(PosF[2]))]
    I = 0.
    J = 0.
    K = 0.
    for i in range(1,len(x_coords)):
        if p_round[0] == x_coords[i]:
            I = i
    for j in range(len(y_coords)):
        if p_round[1] == y_coords[j]:
            J = j
    for k in range(len(z_coords)):
        if p_round[2] == z_coords[k]:
            K = k
                  
    #find velocity corresponding to p_rounded
    Vel = N.array([VelF[0][K][J][I], VelF[1][K][J][I], VelF[2][K][J][I]])
    
    return Vel
    
def findPosIjk(p,PosF):
    """
    function that looks for the indices of the rounded position of a point
    ---Inputs---
    p - point where the indices will be approximated
    PosF - position field where 'p' is located in. [X,Y,Z]
    ---Outputs---
    ijk - the rounded indices of p
    ---
    C.Losada de la Lastra 2015
    """
    #round coordinates of p
    p_round = [round(p[i]) for i in range(len(p))]
    
    #find indices corresponding to p_rounded
    x_coords = PosF[0][0][0]
    y_coords = [PosF[1][0][i][0] for i in range(len(PosF[1][0]))]
    z_coords = [PosF[2][i][0][0] for i in range(len(PosF[2]))]
    I = 0.
    J = 0.
    K = 0.
    for i in range(len(x_coords)):
        if p_round[0] == x_coords[i]:
            I = i
    for j in range(len(y_coords)):
        if p_round[1] == y_coords[j]:
            J = j
    for k in range(len(z_coords)):
        if p_round[2] == z_coords[k]:
            K = k
    ijk = [I,J,K]
    return ijk