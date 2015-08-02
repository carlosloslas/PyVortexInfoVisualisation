# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 10:43:53 2014
@author: closadalastra
"""
"""
Script that joins the generation of the vel. field and the visualisation tools
vel. field: VTK file
visual. tools: lambda2, seeding tecnique, geometricVortex
"""
import numpy as N

#functions to read the field from the vtk image
from VTK_ReadVf3d import VTK_ReadVf3d

#functions to visualise the field
from lambda2 import lambda2_3d
from geometricVortex import geoVortex

##filename of, and format of, VTK file where the velocity field is contained.
velocityfile = 'vrtx04_24_15_12rad4core60pnts.vtu'
    #'vr04_17_15_7rad2core50pnts.vtu' slightly bigger vortex
    #'vr04_17_15_6rad2core60pnts.vtu' Small single vortex with core
    #'vrtx04_24_15_12rad4core60pnts.vtu' Best singe vortex with core
    #'vrtx02_26_15_pair3core20r52pnts.vtu'
    #'vrtx02_26_15_2coreWake52pnts.vtu'
    #'vortex02_18_15_pair52pnts.vtu' , vortex pair.
    #'vortex5wake04_22_15_6rad2core60pnts.vtu' , vortex wake.
    #'vortex5wake04_20_15_6rad2core60pnts.vtu'
vtkformat = 'vtu'
###---read VTK---###
data = VTK_ReadVf3d(velocityfile,vtkformat)
[X,Y,Z,U,V,W] = data
PosF = [X,Y,Z]
VelF = [U,V,W]
#xmax = N.amax(X)
#xmin = N.amin(X)
#ymax = N.amax(Y)
#ymin = N.amin(Y)
#zmax = N.amax(Z)
#zmin = N.amin(Z)
print('Data read from '+velocityfile)


###---Calcs---###

###lambda2_3d(X,Y,Z,U,V,W):
L2 = lambda2_3d(X,Y,Z,U,V,W)


#geoVortex(L2,PosF,VelF)
CircF = geoVortex(L2,PosF,VelF)


print 'Visualisation ended'
