# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 10:43:53 2014
@author: closadalastra
"""
"""
Script that joins the generation of the vel. field and the visualisation tools
velocity field: VTK file
visualisation tools: lambda2, seeding tecnique, geometricVortex
"""
import numpy as N

#functions to read the field from the vtk image
from VTK_ReadVf3d import VTK_ReadVf3d

#functions to visualise the field
from mayavi import mlab
from lambda2 import lambda2_3d
from geometricVortex import geoVortex


##filename of, and format of, VTK file where the velocity field is contained.
velocityfile = 'single vortex.vtu'
    #'single vortex.vtu'
    #'vortex wake.vtu'
vtkformat = 'vtu'


###---read VTK---###
data = VTK_ReadVf3d(velocityfile,vtkformat)
[X,Y,Z,U,V,W] = data
PosF = [X,Y,Z]
VelF = [U,V,W]

print('=> Data read from '+velocityfile)


###---Calcs---###

#Lambda 2 definition of a vortex core.
###lambda2_3d(X,Y,Z,U,V,W):
L2 = lambda2_3d(X,Y,Z,U,V,W)
print('=> Lambda_2 field calculated')

#Geometric and vortex property calculations.
#geoVortex(L2,PosF,VelF):
CircF = geoVortex(L2,PosF,VelF)
print '=> Visualisation ended'
