# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:05:15 2015

@author: closadalastra
"""
import numpy as N
import math as M
from numpy import linalg as LA
from mayavi import mlab

#@mlab.show
def geoVortex(L2,PosF,VelF):
    """
    function that takes the L2 field, seeds it, and calculates an aproximate circulation
    ---Inputs---
    L2 - the lambda2 array
    L2Vec - the lambda2 eigenvector array
    PosF - Position field
    ---Output---
    
    ---
    C.Losada de la Lastra 2015
    """
    from lambda2 import L2min, cropL2min
    from seeding import seed_time
    from circulationVortex import circulationPoint,circulationFieldUpdate,hasPointCirculation
    
    X,Y,Z = PosF 
    
    fig1 = mlab.figure(1, size=(1000, 500), fgcolor=(1, 1, 1),\
        bgcolor=(0.5, 0.5, 0.5))     
    #Lambda2Surf = mlab.contour3d(Z,Y,X,L2,figure=fig1,contours=3,transparent=True)
#    Lambda2Legend = mlab.scalarbar(object=Lambda2Surf, title='Lambda_2',\
#        orientation='vertical', nb_labels=None, nb_colors=None, label_fmt=None)
    print('displayed L2 isosurface')
    
    ##---1---## L2 min points
    #L2min(PosF,L2)
    dataL2min = L2min(PosF,L2)
    
    ##---2---## Crop L2min points
    #cropL2min(dataL2min,PosF,method='negData')
    #L2minData_crop = cropL2min(dataL2min,PosF,method='auto')
    L2minData_crop = cropL2min(dataL2min,PosF,method='negData')
    print('obtained and croped L2min data')
    seedsL2min_l2 = L2minData_crop[0]
    seedsL2min_p = L2minData_crop[1]
    seedsL2min_ijk = L2minData_crop[2]
    
    ##---3---## Eval in time L2min seeds 
    seedsL2min_paths_time = []
    seedsL2min_vels_time = []
    seedsL2min_circ = []
    seedsL2min_swirlVects = []
    for i in range(len(seedsL2min_p)):
        p0 = seedsL2min_p[i]
        ijk = seedsL2min_ijk[i]
        t = 4
        seed_pi2 = 0
        while seed_pi2==0:
            #seed_time(p0,ijk,PosF,VelF,t)
            path_s,vel_s = seed_time(p0,ijk,PosF,VelF,t)
            #findProbeAng(path)
            #crit2pi(thetas)
            pi2 = crit2pi(findProbeAng(path_s))
            if pi2 == True:
                seed_pi2 = 1
                seedsL2min_paths_time.append(path_s)                
                seedsL2min_vels_time.append(vel_s)
                #circulationPoint(vel_seed_time)
                seedsL2min_circ.append(circulationPoint(vel_s))
                #averageNormal(path)
                seedsL2min_swirlVects.append(averageNormal(path_s))
                
                #plot in MAYAVI
                mlab.plot3d(path_s[2],path_s[1],path_s[0],figure=fig1)                
                
            else:
                
                t+=1
    ##--4--## Circulation of seedsL2min into circulation field
    CircF0 = N.zeros_like(X,dtype=float)
    #circulationFieldUpdate(pnts,pnts_ijk,pnts_circ,pnts_vect_swirl,CircF,PosF,L2)
    pnts = seedsL2min_p
    pnts_ijk = seedsL2min_ijk
    pnts_circ = seedsL2min_circ
    pnts_vect_swirl = seedsL2min_swirlVects
    CircF = circulationFieldUpdate(pnts,pnts_ijk,pnts_circ,pnts_vect_swirl,CircF0,PosF,L2)
#    Vmin = 
    CircField = mlab.contour3d(Z,Y,X,CircF,contours=11,vmin=min(CircF.flatten()),vmax=-min(CircF.flatten()))
    CircFieldLegend = mlab.scalarbar(object=CircField, title='Circulation',\
        orientation='vertical', nb_labels=5, nb_colors=None, label_fmt=None)
    
    return CircF  

  
    ##--5--## Integrate L2 surf from each seedL2min till nextL2point has a circulation
##    for i in range(0,1):#len(seedsL2min_p)):
#    i=-1
#    print 'starting seed '+str(i)+' out of '+str(len(seedsL2min_p))
#    p0 = seedsL2min_p[i]
#    p0_ijk = seedsL2min_ijk[i]
#    #p0_path = seedsL2min_paths_time[i]
#    p0_vect = seedsL2min_swirlVects[i]
#    p0_circ = seedsL2min_circ[i]
#    
#    L2points = []
#    L2points.append(p0)
#    L2points_vect = []
#    L2points_vect.append(p0_vect)
#    L2points_ijk = []
#    L2points_ijk.append(p0_ijk)
#    
#    L2points_circ = []
#    L2points_circ.append(p0_circ)
#    
#    next_L2point_circ = 0
#    n =  0
#    t = 8
#    dist = 20 #-> distance to next L2point <-#
#    while next_L2point_circ < 4:
#        L2point0 = L2points[n]
#        L2point_ijk0 = L2points_ijk[n]
#        L2point_vect0 = L2points_vect[n]
#        #findNextL2point_3rdAprox(p0,ijk,dist,vect_dir,L2,PosF,VelF)
#        p1,ijk1 = findNextL2point_3rdAprox(L2point0, L2point_ijk0,\
#                        dist, L2point_vect0, L2,PosF,VelF)
#        #hasPointCirculation(CircF,ijk_p)
#        if hasPointCirculation(CircF,ijk1) == True:
#            next_L2point_circ += 1
#        else:
#            pass
#        p1_pi2=0
#        while p1_pi2 == 0:
#            #seed_time(p0,ijk,PosF,VelF,t)
#            path_p1,vel_p1 = seed_time(p1,ijk1,PosF,VelF,t)
#            #findProbeAng(path)
#            #crit2pi(thetas)
#            pi2 = crit2pi(findProbeAng(path_s))
#            if pi2 == True:
#                p1_pi2=1
#                L2points.append(p1)
#                L2points_ijk.append(ijk1)
#                #averageNormal(path)
#                L2points_vect.append(averageNormal(path_p1))
#                #circulationPoint(vel_seed_time)
#                L2points_circ.append(circulationPoint(vel_p1))
#                n+=1
#                    
#                #plot in MAYAVI
#                mlab.plot3d(path_p1[2],path_p1[1],path_p1[0],figure=fig1)                        
#                    
#            else:
#                t+=2
##        #add calculated circulation to circulation field
##        if len(L2points) < 4:
##            pass
##        else:
##            #circulationFieldUpdate(pnts,pnts_ijk,pnts_circ,pnts_vect_swirl,CircF,PosF,L2)
##            CircF = circulationFieldUpdate(L2points,L2points_ijk,L2points_circ,\
##                                            L2points_vect,CircF,PosF,L2)
##        print 'Done seed '+str(i)+' out of '+str(len(seedsL2min_p))
#        
#    CircField = mlab.contour3d(Z,Y,X,CircF,contours=10)
#    CircFieldLegend = mlab.scalarbar(object=CircField, title='Circulation',\
#        orientation='vertical', nb_labels=None, nb_colors=None, label_fmt=None)
#        
#    ###---DONE?---###
#    mlab.show()
#    return CircF

    #CircF = N.zeros_like(X,dtype=float)
    
#    seeds_l2 = L2minData_crop[0]
#    seeds_p = L2minData_crop[1]
#    seeds_ijk = L2minData_crop[2]
#    
#    seeds_test_l2 = []
#    seeds_test_p = []
#    seeds_test_ijk = []
#    
#    seeds_test_l2.append(seeds_l2[-1])
#    seeds_test_p.append(seeds_p[-1])
#    seeds_test_ijk.append(seeds_ijk[-1])
#    
#    seeds_time_test = []
#    seeds_time_test_mayavi = []
#    seeds_new_mayavi = []
#    
#    seeds_test_circulation = []
#    seeds_test_circulation.append()
#    seeds_test_swirl = []
    
#    ###eval first point
#    i=0
#    t=8
#    n_s=0
#    path_seed_time,vel_seed_time = seed_time(seeds_test_p[i],\
#                                                seeds_test_ijk[i],PosF,VelF,t)
#    seeds_time_test.append(path_seed_time)
#    
#    #findProbeAng(path)
#    probeAng = findProbeAng(path_seed_time)
#    #crit2pi(thetas)
#    pi2 = crit2pi(probeAng)
#    if pi2 == True:
#        #plot the path
#        sx,sy,sz = path_seed_time
#        seeds_time_test_mayavi.append(mlab.plot3d(sz,sy,sx,figure=fig1))
##            seeds_time_test_mayavi.append(mlab.plot3d(sz,sy,sx))
##            mlab.draw(figure=fig1)
#        #mlab.show(fig1)
#        
#        #find next L2 point
#        dist = 20
#        #averageNormal(path)
#        vect_dir = averageNormal(path_seed_time)
#        #findNextL2point_3rdAprox(p0,ijk,dist,vect_dir,L2,PosF,VelF):
#        [p1,ijk1] = findNextL2point_3rdAprox(seeds_test_p[i],seeds_test_ijk[i],\
#                                    dist,vect_dir,L2,PosF,VelF)
#        #seeds_new_mayavi.append(mlab.points3d((p1[2],p1[1],p1[0]),figure=fig1,scale_factor=0.2))
##            seeds_new_mayavi.append(mlab.points3d((p1[2],p1[1],p1[0]),scale_factor=0.2))
#        #origin = mlab.points3d(0,0,0,scale_factor=0.1)
#        #mlab.draw(fig1)
#       # mlab.show(fig1)
#        seeds_test_p.append(p1)
#        seeds_test_ijk.append(ijk1)
#        i+=1
#        n_s += 1
#        t=8
#        
#    else:
#        t +=2    
#    pnts=[seeds_test_p[0]]
#    pnts_ijk=[seeds_test_ijk[0]]
#    #circulationPoint(vel_seed_time)
#    pnts_circ=[circulationPoint(path_seed_time)]
#    pnts_vect_swirl = [vect_dir]
#    CircF1 = circulationFieldUpdate(pnts,pnts_ijk,pnts_circ,pnts_vect_swirl,CircF,PosF,L2)
#    
#    CircField = mlab.contour3d(Z,Y,X,CircF1,contours=10)
#    CircFieldLegend = mlab.scalarbar(object=CircField, title='Circulation',\
#        orientation='vertical', nb_labels=None, nb_colors=None, label_fmt=None)
    
##    geo = 0
#    t = 8
#    i = 0
##    n_s=0
#    #n_s_total = 60  #number of seeeds to calc geometric stuff
#    point_circ = 0    
#    print('Starting geometric calculations...')
##    l2pnts=0
##    while l2pnts < len()
#    while point_circ < 4:
#        #calc seed in time
#        #seed_time(p0,ijk,PosF,VelF,t)
#        path_seed_time,vel_seed_time = seed_time(seeds_test_p[i],\
#                                                seeds_test_ijk[i],PosF,VelF,t)
#        seeds_time_test.append(path_seed_time)
#        
#        #findProbeAng(path)
#        probeAng = findProbeAng(path_seed_time)
#        #crit2pi(thetas)
#        pi2 = crit2pi(probeAng)
#        if pi2 == True:
#            #plot the path
#            sx,sy,sz = path_seed_time
#            seeds_time_test_mayavi.append(mlab.plot3d(sz,sy,sx,figure=fig1))
#            
#            seeds_test_circulation.append(circulationPoint(vel_seed_time))
#            
#            #find next L2 point
#            dist = 20
#            #averageNormal(path)
##            vect_dir = averageNormal(path_seed_time)
#            seeds_test_swirl.append(averageNormal(path_seed_time))
#            #findNextL2point_3rdAprox(p0,ijk,dist,vect_dir,L2,PosF,VelF):
#            [p1,ijk1] = findNextL2point_3rdAprox(seeds_test_p[i],seeds_test_ijk[i],\
#                                        dist,seeds_test_swirl[i],L2,PosF,VelF)
#            if hasPointCirculation(CircF,ijk1) == True:
#                point_circ += 1
#            else:
#                pass
#            
#            #seeds_new_mayavi.append(mlab.points3d((p1[2],p1[1],p1[0]),figure=fig1,scale_factor=0.2))
##            seeds_new_mayavi.append(mlab.points3d((p1[2],p1[1],p1[0]),scale_factor=0.2))
#            #origin = mlab.points3d(0,0,0,scale_factor=0.1)
#            #mlab.draw(fig1)
#           # mlab.show(fig1)
#            seeds_test_p.append(p1)
#            seeds_test_ijk.append(ijk1)
#            i+=1
#            #n_s += 1
#            t=5
#                        
#        else:
#            t +=2
#    CircF = circulationFieldUpdate(seeds_test_p,seeds_test_ijk,seeds_test_circulation,seeds_test_swirl,CircF,PosF,L2)
#    
#    CircField = mlab.contour3d(Z,Y,X,CircF,contours=10)
#    CircFieldLegend = mlab.scalarbar(object=CircField, title='Circulation',\
#        orientation='vertical', nb_labels=None, nb_colors=None, label_fmt=None)
#    
#    mlab.show()

def findNextL2point_4thAprox(p0,ijk,dist,vect_dir,L2,PosF,VelF):
    """
    function that finds the next point in the field to carry out 
    the geometric calculation using a fourth order approximation.
    ---Inputs---
    p0 - previous point in 'L2'
    ijk - indices of 'p' in 'L2' and 'PosF'
    dist - distance from 'p' to the new point
    vect_dir - direction in which to find the next point
    L2 - Lambda2 scalar field. 3d N.array
    PosF - position field corresponding to 'L2' [X,Y,Z]
    VelF - velocity field corresponding to 'L2' [U,V,W]
    ---Output---
    point - the next point in 'PosF' [[xp,yp,zp],[ip,jp,kp]]
    ---
    C.Losada de la Lastra 2015
    """
    from seeding import seed_time, findPosIjk
    
    #findNextL2point(p0,ijk,dist,vect_dir,L2,PosF):
    p1,ijk1 = findNextL2point(p0,ijk,dist,vect_dir,L2,PosF)
    
    
    #second order approx
    p1_2pi = 0
    t=5
    while p1_2pi == 0:
        #seed_time(p0,ijk,PosF,VelF,t): 
        p1_time,p1_vel = seed_time(p0,ijk,PosF,VelF,t)
        #findProbeAng(path)
        #crit2pi(thetas)
        pi2 = crit2pi(findProbeAng(p1_time))
        if pi2 == True:
            p1_2pi = 1
            #averageNormal(path)
            p1_norm_vect = averageNormal(p1_time)
        else:
            t+=2
                
    dist_p0p1 = LA.norm(p1-p0)
    p2 = p0 + (0.5*dist_p0p1*vect_dir) + (0.5*dist_p0p1*p1_norm_vect)
    #findPosIjk(p,PosF):
    ijk2 = i2,j2,k2 = findPosIjk(p2,PosF)
    
    #third order aprox
    p2_2pi = 0
    t=5
    while p2_2pi == 0:
        #seed_time(p0,ijk,PosF,VelF,t): 
        p2_time,p2_vel = seed_time(p2,ijk2,PosF,VelF,t)
        #findProbeAng(path)
        #crit2pi(thetas)
        pi2 = crit2pi(findProbeAng(p2_time))
        if pi2 == True:
            p2_2pi = 1
            #averageNormal(path)
            p2_norm_vect = averageNormal(p2_time)
        else:
            t+=2   
    
    dist_p0p1 = LA.norm(p1-p0)
    p3 = p0 + ((1/3.)*dist_p0p1*vect_dir) + ((1/3.)*dist_p0p1*p1_norm_vect) +\
            ((1/3.)*dist_p0p1*p2_norm_vect)
    #findPosIjk(p,PosF):
    ijk3 = i3,j3,k3 = findPosIjk(p3,PosF)

    #fourth order aprox
    p3_2pi = 0
    t=5
    while p3_2pi == 0:
        #seed_time(p0,ijk,PosF,VelF,t): 
        p3_time,p3_vel = seed_time(p3,ijk3,PosF,VelF,t)
        #findProbeAng(path)
        #crit2pi(thetas)
        pi2 = crit2pi(findProbeAng(p3_time))
        if pi2 == True:
            p2_2pi = 1
            #averageNormal(path)
            p3_norm_vect = averageNormal(p3_time)
        else:
            t+=2       
    
    dist_p0p1 = LA.norm(p1-p0)
    p4 = p0 + (0.25*dist_p0p1*vect_dir) + (0.25*dist_p0p1*p1_norm_vect) +\
            (0.25*dist_p0p1*p2_norm_vect) + (0.25*dist_p0p1*p3_norm_vect)
    #findPosIjk(p,PosF):
    i4,j4,k4 = findPosIjk(p4,PosF)
    return [p4,[i4,j4,k4]]

def findNextL2point_3rdAprox(p0,ijk,dist,vect_dir,L2,PosF,VelF):
    """
    function that finds the next point in the field to carry out 
    the geometric calculation using a third order approximation.
    ---Inputs---
    p0 - previous point in 'L2'
    ijk - indices of 'p' in 'L2' and 'PosF'
    dist - distance from 'p' to the new point
    vect_dir - direction in which to find the next point
    L2 - Lambda2 scalar field. 3d N.array
    PosF - position field corresponding to 'L2' [X,Y,Z]
    VelF - velocity field corresponding to 'L2' [U,V,W]
    ---Output---
    point - the next point in 'PosF' [[xp,yp,zp],[ip,jp,kp]]
    ---
    C.Losada de la Lastra 2015
    """
    from seeding import seed_time, findPosIjk
    
    #findNextL2point(p0,ijk,dist,vect_dir,L2,PosF):
    p1,ijk1 = findNextL2point(p0,ijk,dist,vect_dir,L2,PosF)
    #second order approx
    p1_2pi = 0
    t=5
    while p1_2pi == 0:
        #seed_time(p0,ijk,PosF,VelF,t): 
        p1_time,p1_vel = seed_time(p0,ijk,PosF,VelF,t)
        #findProbeAng(path)
        #crit2pi(thetas)
        pi2 = crit2pi(findProbeAng(p1_time))
        if pi2 == True:
            p1_2pi = 1
            #averageNormal(path)
            p1_norm_vect = averageNormal(p1_time)
        else:
            t+=2
                
    dist_p0p1 = LA.norm(p1-p0)
    p2 = p0 + (0.5*dist_p0p1*vect_dir) + (0.5*dist_p0p1*p1_norm_vect)
    #findPosIjk(p,PosF):
    ijk2 = i2,j2,k2 = findPosIjk(p2,PosF)
    
    #third order aprox
    p2_2pi = 0
    t=5
    while p2_2pi == 0:
        #seed_time(p0,ijk,PosF,VelF,t): 
        p2_time,p2_vel = seed_time(p2,ijk2,PosF,VelF,t)
        #findProbeAng(path)
        #crit2pi(thetas)
        pi2 = crit2pi(findProbeAng(p2_time))
        if pi2 == True:
            p2_2pi = 1
            #averageNormal(path)
            p2_norm_vect = averageNormal(p2_time)
        else:
            t+=2   
    
    dist_p0p1 = LA.norm(p1-p0)
    p3 = p0 + ((1/3.)*dist_p0p1*vect_dir) + ((1/3.)*dist_p0p1*p1_norm_vect) +\
            ((1/3.)*dist_p0p1*p2_norm_vect)
    #findPosIjk(p,PosF):
    i3,j3,k3 = findPosIjk(p3,PosF)
    return [p3,[i3,j3,k3]]
  
def findNextL2point_2ndAprox(p0,ijk,dist,vect_dir,L2,PosF,VelF):
    """
    function that finds the next point in the field to carry out 
    the geometric calculation using a second order approximation.
    ---Inputs---
    p0 - previous point in 'L2'
    ijk - indices of 'p' in 'L2' and 'PosF'
    dist - distance from 'p' to the new point
    vect_dir - direction in which to find the next point
    L2 - Lambda2 scalar field. 3d N.array
    PosF - position field corresponding to 'L2' [X,Y,Z]
    VelF - velocity field corresponding to 'L2' [U,V,W]
    ---Output---
    point - the next point in 'PosF' [[xp,yp,zp],[ip,jp,kp]]
    ---
    C.Losada de la Lastra 2015
    """
    from seeding import seed_time, findPosIjk
    
    #findNextL2point(p0,ijk,dist,vect_dir,L2,PosF):
    p1,ijk1 = findNextL2point(p0,ijk,dist,vect_dir,L2,PosF)
    p1_2pi = 0
    t=5
    while p1_2pi == 0:
        #seed_time(p0,ijk,PosF,VelF,t): 
        p1_time,p1_vel = seed_time(p0,ijk,PosF,VelF,t)
        #findProbeAng(path)
        #crit2pi(thetas)
        pi2 = crit2pi(findProbeAng(p1_time))
        if pi2 == True:
            p1_2pi = 1
            #averageNormal(path)
            p1_norm_vect = averageNormal(p1_time)
        else:
            t+=2
                
    dist_p0p1 = LA.norm(p1-p0)
    p2 = p0 + (0.5*dist_p0p1*vect_dir) + (0.5*dist_p0p1*p1_norm_vect)
    #findPosIjk(p,PosF):
    i2,j2,k2 = findPosIjk(p2,PosF)
    return [p2,[i2,j2,k2]]
    
def findNextL2point(p0,ijk,dist,vect_dir,L2,PosF):
    """
    function that finds the next point in the field to carry out 
    the geometric calculation
    ---Inputs---
    p0 - previous point in 'L2'
    ijk - indices of 'p' in 'L2' and 'PosF'
    dist - distance from 'p' to the new point
    vect_dir - direction in which to find the next point
    L2 - Lambda2 scalar field. 3d N.array
    PosF - position field corresponding to 'L2' [X,Y,Z]
    ---Output---
    point - the next point in 'PosF' [[xp,yp,zp],[ip,jp,kp]]
    ---
    C.Losada de la Lastra 2015
    """
    from seeding import findPosIjk
    
    i0,j0,k0 = ijk
    l2_p0 = L2[k0][j0][i0]
    
    pnt_found = 0
    while  pnt_found == 0:
        p1 = p0+vect_dir * dist
        p1_r = N.array([round(p1[i]) for i in range(len(p1))])
        #findPosIjk(p,PosF):
        i1,j1,k1 = findPosIjk(p1_r,PosF)
        l2_p1 = L2[k1][j1][i1]
        if l2_p1 < 0:            
            if float(truncate(l2_p0,1))  == round(l2_p1,1):
                pnt_found = 1
            else:
                dist *= 0.9
        else:
            dist *= 0.8
            
    return [p1,[i1,j1,k1]]
    
def averageNormal(path):
    """
    function that calculates the average direction in which to find the 
    next seed in order to continue applying the geometirc calculations 
    ---Inputs---
    path - path of the previous seed.
    ---Output---
    vect_dir - vector in the direction in which to find the next seed
    ---
    C.Losada de la Lastra 2015
    """    
    vects = N.array([N.array([path[0][i],path[1][i],path[2][i]])-N.array([\
                    path[0][i-1],path[1][i-1],path[2][i-1]]) for i in range(\
                    1,len(path[0]))])
    n_vects = N.array([N.cross(vects[i-1],vects[i]) for i in range(1,len(vects))])
    len_n_vects = 0.0
    for i in range(len(n_vects)):
        if (n_vects[i] == N.array([0.,0.,0.])).all():
            pass
        else:
            len_n_vects += 1
    
    vect_dir = sum(n_vects)/len_n_vects
    
    return vect_dir
    
def crit2pi(thetas):
    """
    function that confirms weather or not the 2pi criterion has been met.
    ---Inputs---
    thetas - probe angles between the points as produced by the findAng function.
    ---Output---
    crit - logical true or false, depending on weather the criterion has been met.
    ---
    C.Losada de la Lastra 2015
    """
    angle = 0 
    for i in range(len(thetas)):
        angle += thetas[i]
        if angle >= 2*M.pi:
            return True
        else:
            pass
    return False

def findProbeAng(path):
    """
    function that calculates the angle between the points allong the 
    path of a seed
    ---Inputs---
    path - points that form the path of a seed with time 
    [[x0,y0,z0],[x1,y1,z1],...,[xt,yt,zt]]
    ---Output---
    thetas - angle allong the path of the seed.
    
    ---
    C.Losada de la Lastra 2015
    """
    thetas = N.zeros([len(path[0])-2])
    
    for i in range(1,len(path[0])-1):
        p_0 = N.array([path[0][i-1],path[1][i-1],path[2][i-1]])
        p_1 = N.array([path[0][i],path[1][i],path[2][i]])
        p_2 = N.array([path[0][i+1],path[1][i+1],path[2][i+1]])
        thetas[i-1] = findAng(p_0,p_1,p_2) 
    
    return thetas

def findAng(p0,p1,p2):
    """
    function that calculates the angle between the two vectors formed by p0p1 and p1p2
    ---Inputs---
    p0,p1,p2 - three points allong the trajectory of the seed
    ---Output---
    theta - angle between the two vectors IN RADIANS.
    ---
    C.Losada de la Lastra 2015
    """
    v0 = p1 - p0
    v1 = p2 - p1
    v0_u = v0/LA.norm(v0)
    v1_u = v1/LA.norm(v1)
    
    vu_dot = N.dot(v0_u,v1_u)    
    try:
        theta = M.acos(vu_dot) #IN RADIANS
    except ValueError:
        print 'cos of angle between vectors is out of range...rounding...'
        vu_dot_r = round(vu_dot,14)
        return M.acos(vu_dot_r)
    else:            
        return theta
        
def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])
    