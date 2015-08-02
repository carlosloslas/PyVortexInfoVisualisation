# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 15:46:11 2015

@author: closadalastra
"""
import numpy as N
import vtk
from vtk.util.numpy_support import vtk_to_numpy as vtk2np

def VTK_ReadVf3d(name,VTKformat):
    """
    VTK_ReadVf3d(name,VTKformat):
    function that calls the appropiate developed vtk reader for the appropiate 
    vtk format.
    outputs a dictionary with the arrays of the points, scalars, and vectors.
    numpy.array([points,scalars,vectors])
    -name => name of the vtk file to read. 
    -VTKformat => Format of the file to be read.
        currently available formats:
            vtu: binary unstructured grid file.
            vti: to be developed
            vts: to be developed
    """
    #import vf3d_vtu
    
    if VTKformat == 'vtu':
        res_dict = vtu_readvf3d(name)
    elif VTKformat == None:
        print 'Please select a VTK format'
    else:
        print 'The selected format has not been developed yet'
    return res_dict

#---.vtu format---#
def vtu_readvf3d(name):
        
    #generate reader and read file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name)
    reader.Update()
    vtkdata = reader.GetOutput()
    
    #Obtain vectors and points from data in file
    nPoints = vtkdata.GetNumberOfPoints()
    pointData = vtkdata.GetPointData()
    vectorData = pointData.GetVectors()
    scalarData = pointData.GetScalars()    
    try:
        vect_list = vtk2np(vectorData)
    except AttributeError:
        print 'there is no vector data'
        vect_list = None
    try:
        scal_list = vtk2np(scalarData)
    except AttributeError:
        #print 'there is no scalar data'
        scal_list = None
    pnts_list = N.array([vtkdata.GetPoint(i) for i in range(nPoints)])

    #rearange data to 3d matices    
    bounds = vtkdata.GetBounds()#(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax) 
    [X,Y,Z,U,V,W] = ijk(pnts_list,bounds,vect_list) 
    
    return [X,Y,Z,U,V,W]
    
 #---.vti format---# to be developed

 #---.vts format---# to be developed

###---util functions---###
#1.ikj(points,bounds,vectors)
#2.dataVtk_3dMatrix(points,bounds,vectors)
#
def ijk(points,bounds,vectors):
    """
    turns the points and vectors to 3d arrays of X Y Z U V W.
    returns an array with the 3d arrays associated to X Y Z U V W.
    -points => list of the coordinates of the poitns where the data is located.
    -bounds => bounds of the data.(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
    -vectors => vector data of the field at the 'points'
    """
    #this functions assumes EQUAL length axis of +ve and -ve directions.
    #asign variables
    (xmin,xmax,ymin,ymax,zmin,zmax) = bounds
    #generate the output arrays
    grid3d = N.mgrid[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1]
    U = N.zeros_like(grid3d[0])
    V = N.zeros_like(grid3d[0])
    W = N.zeros_like(grid3d[0])
    X = N.zeros_like(grid3d[0])
    Y = N.zeros_like(grid3d[0])
    Z = N.zeros_like(grid3d[0])
    #loop and rearange
    for n in range(len(points)):
        i = points[n][0]+xmax
        j = points[n][1]+ymax
        k = points[n][2]+zmax
        X[k][j][i] = points[n][0]
        Y[k][j][i] = points[n][1]
        Z[k][j][i] = points[n][2]
        U[k][j][i] = vectors[n][0]
        V[k][j][i] = vectors[n][1]
        W[k][j][i] = vectors[n][2]
        
        float('nan')
        if X[k][j][i] == 'nan':
            print k,j,i
        elif Y[k][j][i] == 'nan':
            print k,j,i
        elif Z[k][j][i] == 'nan':
            print k,j,i
        elif U[k][j][i] == 'nan':
            print k,j,i
        elif V[k][j][i] == 'nan':
            print k,j,i
        elif W[k][j][i] == 'nan':
            print k,j,i
        
    return [X,Y,Z,U,V,W]

def dataVtk_3dMatrix(points,bounds,vectors):
    """
    Function that turns a vtk output formated data to 3d field matrix data
    from [(x1,y1,z1),...,(xn,yn,zn)]
    to [[[[x1,y1,z1],[...],[x3,y1,z1]],[[x1,y2,z1],[...],[...]],[[x1,y3,z1],[...],[...]]]
    ,[[[x1,y1,z2],[...],[...]],[...],[...]] , [.........]]    
    -points => list of the coordinates of the poitns where the data is located.
    -bounds => bounds of the data.(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
    -vectors => vector data of the field at the 'points'
    """
    #asign variables
    (xmin,xmax,ymin,ymax,zmin,zmax) = bounds
        
    #generate the output arrays
    grid3d = N.mgrid[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1]
    pnts3d = N.zeros_like(grid3d[0],dtype= N.ndarray)
    vect3d = N.zeros_like(grid3d[0],dtype= N.ndarray)
    
    #loop and rearange
    for i in range(len(points)):
        x_t = points[i][0]
        y_t = points[i][1]
        z_t = points[i][2]
        pnts3d[z_t+zmax][y_t+ymax][x_t+xmax] = points[i]
        vect3d[z_t+zmax][y_t+ymax][x_t+xmax] = vectors[i]
        
    return {'points':pnts3d,'vectors':vect3d}


#VTK_ReadVf3d('vortex02_4_15.vtu','vtu')
#VTK_ReadVf3d('vtktest.vtu','vtu')
