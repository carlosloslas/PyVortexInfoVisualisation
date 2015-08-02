# -*- coding: utf-8 -*-
"""
FEEG3003 Individual Project,
3D visualisation of cuantitative vortex information.
"""

def vField3d(C, r, Nvec, np, Rcore=None):
    """
    Function that generates the 3D velocity field induced by a nunmber
    of 3D vortex rings.
    ---Inputs---
    C - centroid of vortex ring/s
    r - radious of ring/s
    Nvec - normal vector to the ring/s
    np - number of points that the ring is defined by
    Rcore - radious of the core of the vortex, in number of units.
    ---
    C.Losada de la Lastra 2014
    """
    import numpy as N    
    from mayavi import mlab
    
    from circumference import circumference
    from tangent import tangent

    from vortexLine import vortexLine3  
    from vortexCore import vortexCore3
    from vF3d_VTK import vF3d_VTK
    
    gamma=1 #strength of vortex rings
        
    #size of the field in relation to size of rings
    x_coords= N.array([C[i+1][0] for i in range(C[0][0])])
    y_coords= N.array([C[i+1][1] for i in range(C[0][0])])
    z_coords= N.array([C[i+1][2] for i in range(C[0][0])])
        
    f=1.25 #parameter controling the size of the field    
    xmin= min(x_coords)- f*r
    xmax= max(x_coords)+ f*r
    ymin= min(y_coords)- f*r
    ymax= max(y_coords)+ f*r
    zmin= min(z_coords)- f*r
    zmax= max(z_coords)+ f*r
    
    #generation of arrays where the coordinates of each point arelocated
                                #and velocity components will be stored
    #array of x,y,z coordinates
    [Z,Y,X] = N.mgrid[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1]

    #Points constituting the Vortex Ring/s
    #circumference(c,r,n,np)
    C_out=N.array([circumference(C[i+1],r,Nvec[i+1],np) for i in range(C[0][0])])
    
    #VR points which constitute the vortex rings. 
    #   points with which the tangent vectors are calculated (VR = VortexRing)
    #VR_iVel points where the induced velocity is calculated about(VortexRing_inducedVelocity)
    VR = N.array([C_out[i][0] for i in range(C[0][0])])
    VR_iVel = N.array([C_out[i][1] for i in range(C[0][0])])
    
    #Vectors Tangent to the Ring/s
    #tangent(C)
    TVR = N.array([tangent(VR[i]) for i in range(len(VR))])
    
    #Generating Velicity Field
    if Rcore == None:
        #vortexLine3(PosF,C,TVR,VR_iVel,gamma):
        [U,V,W] = vortexLine3([X,Y,Z],C,TVR,VR_iVel,gamma)
    else:
        #vortexCore3(PosF,C,TVR,VR_iVel,Rcore,gamma):
        [U,V,W] = vortexCore3([X,Y,Z],C,TVR,VR_iVel,Rcore,gamma)
    
    ###---VISUALISE IN MAYAVI---###
    vis = 1 # ON/OFF
    if vis == 1:
        fig = mlab.figure(1, size=(1000, 500), fgcolor=(1, 1, 1),\
            bgcolor=(0.5, 0.5, 0.5))
        print('Plot generated') 
        #vortex rings
        vr_plt = N.array([mlab.plot3d(VR[i][0],VR[i][1],VR[i][2]) for i in range(C[0][0])]) 
        #velocity Field
        velField =mlab.quiver3d(X,Y,Z,U,V,W,figure=fig,name='VelocityField')
        mlab.show()
    else:
        pass
        
    ###---Write to .vtu file---###
    vtu = 1 # ON/OFF
    if vtu == 1:
        #fname = 'vortex02_26_15_2core32pnts'
        #fname = 'vrtx02_26_15_2coreWake52pnts'    
        #fname = 'vrtx02_26_15_pair3core20r52pnts'
        fname = 'vrtx04_24_15_12rad4core60pnts'    
        #fname = 'vortex5wake04_24_15_6rad2core60pnts'
        fformat = 'vtu'

        vField = N.array([X,Y,Z,U,V,W])                    
        vF3d_VTK(vField,fname,fformat)
        print('vtk image written') 
    else:
        pass
        
    return
    
import numpy as N
vField3d(C=N.array([[1,0,0],[0,0,0]]), r=12, Nvec=N.array([[1,0,0],[1,0,0]]), np=60, Rcore=4)
#vField3d(C=N.array([[1,0,0],[0,0,0]]), r=9, Nvec=N.array([[1,0,0],[1,0,0]]), np=32, Rcore=2)
#vField3d(C=N.array([[5,0,0],[-24,0,0],[-12,0,0],[0,0,0],[12,0,0],[24,0,0]]), r=9, Nvec=N.array([[5,0,0],[1,1,0],[-1,1,0],[1,1,0],[-1,1,0],[1,1,0]]), np=52, Rcore=2)
#vField3d(C=N.array([[2,0,0],[-30,0,0],[30,0,0]]), r=20, Nvec=N.array([[2,0,0],[1,1,0],[-1,1,0]]), np=52, Rcore=3)
#vField3d(C=N.array([[1,0,0],[0,0,0]]), r=3, Nvec=N.array([[1,0,0],[1,0,0]]), np=8)
#vField3d(C=N.array([[1,0,0],[0,0,0]]), r=3, Nvec=N.array([[1,0,0],[1,0,0]]), np=4,Rcore=1)
#vField3d(C=N.array([[1,0,0],[0,0,0]]), r=7, Nvec=N.array([[1,0,0],[1,0,0]]), np=50, Rcore=2)
#vField3d(C=N.array([[1,0,0],[0,0,0]]), r=6, Nvec=N.array([[1,0,0],[1,0,0]]), np=60, Rcore=2)
#vField3d(C=N.array([[5,0,0],[-60,0,0],[-30,0,0],[0,0,0],[30,0,0],[60,0,0]]), r=12, Nvec=N.array([[5,0,0],[1,1,0],[-1,1,0],[1,1,0],[-1,1,0],[1,1,0]]), np=60, Rcore=4)


###---utilities---###
'''
#date obtention for file writing
def date():
    import time
    ## dd/mm/yyyy format
    return (time.strftime("%d/%m/%Y"))
    
#Interactive file and format generation:
    import time
    ## dd/mm/yyyy format
    print (time.strftime("%d/%m/%Y"))
    #fname = 'vortex02_12_15_wake32pnts'
    fname = 'vortex02_12_15_single32pnts'    
    fformat = 'vtu'
    write = raw_input('Write vtk file?(y/n): ')
    if write == 'y':
        namecorrect = raw_input('filename = '+str(fname)+'. Do you want to chage the name of the vtk file?(y/n): ')
        if namecorrect == 'y':
            pass
        else:
            fname = raw_input('Type new name for the vtk file: ') 
        formatcorrect = raw_input('fileformat = '+str(fformat)+'. Do you want to chage the format of the vtk file?(y/n): ')
        if formatcorrect == 'y':
            pass
        #else:
        #   fformat = raw_input('Type new name for the vtk file: ')
        vField = N.array([X,Y,Z,U,V,W])                    
        vF3d_VTK(vField,fname,fformat)
        print('vtk image written')    
'''
        