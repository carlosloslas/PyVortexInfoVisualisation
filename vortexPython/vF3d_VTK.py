# -*- coding: utf-8 -*-
"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
"""
from tvtk.api import tvtk, write_data
import numpy as N

def vF3d_VTK(field,name,VTKformat):
    """
    function that calls the appropiate developed vtk writer to write the 
    3d vector field to the desired format.
    ---Inputs---
    field - Vector field which will be added to the vtk file numpy.array([x,y,z,u,v,w])
    path - path with name of file where the field will be written.
    VTKformat - Desired VTK format from the developed.
        vtu: binary unstructured grid file.
    ---
    C.Losada de la Lastra 2015
    """ 
    if VTKformat == 'vtu':
        vf3d_vtu(field,name)
    elif VTKformat == None:
        print 'Please select a VTK format'
    else:
        print 'The selected format has not been developed yet'
    return #nothing, since functions output the written VTK file

#List of Developed Formats
#---.vtu format---#
def vf3d_vtu(field,name):
    """
    Function that wtrites a .vtu file
    ---Inputs---
    field - Vector field which will be added to the vtk file numpy.array([x,y,z,u,v,w])
    name - Name of the .vtu field
    ---
    C.Losada de la Lastra 2015
    """
    [X,Y,Z,U,V,W] = field #3d velocity field
    
    #achieve the correct format
    Pnts = F3d_2_vtkFromat(N.array([X,Y,Z]))    
    velF = F3d_2_vtkFromat(N.array([U,V,W]))    
    #name the vtu file
    if name == None:
        vtu = 'vf3VTU.vtu'
    else:
        vtu = name + '.vtu'
    
    #Generate and write the .vtu file    
    Ugrid = tvtk.UnstructuredGrid()
    Ugrid.points = Pnts
    Ugrid.point_data.vectors = velF
    Ugrid.point_data.vectors.name = 'velocity'
    
    write_data(Ugrid, vtu)
    
    return vtu
    
 #---.vti format---# to be developed

 #---.vts format---# to be developed

###---util functions---###
#1. vf3d_2_vtkFromat(N.array([Fx,Fy,Fz]))
#

def F3d_2_vtkFromat(F3d):
    """
    Function that turns a 3d field F3d in [Fx,Fy,Fz] format to a 
    [[fx1,fy1,fz1], ... ,[fxn,fyn,fzn]] fomat
    Where FX is the z*y*x array corresponding to the x component of the field 
    and fx1 is the x component of the field at the first point
    ---Inputs---
    f3d - field which will be converted by the function.
    ---
    C.Losada de la Lastra 2015
    """
    #asign variables
    [Fx,Fy,Fz] = F3d
    
    #generate the output array
    F3dVTK = N.array([N.zeros(3) for i in range(len(Fx)*len(Fy[0])*len(Fz[0][0]))])
    
    #loop and rearange
    c=0
    for k in range(len(Fz)):
        for j in range(len(Fz[0])):
            for i in range(len(Fz[0][0])):
                #fariables corresponding with the point
                fxn = Fx[k][j][i]
                fyn = Fy[k][j][i]
                fzn = Fz[k][j][i]
                F3dVTK[c] = N.array([fxn,fyn,fzn])
                #update counter            
                c = c+1
                
    return F3dVTK
    
