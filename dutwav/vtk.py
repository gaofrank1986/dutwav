from dutwav.mesh import Mesh
from pyvtk import *
from dutwav.__mesh_core import POS
from copy import copy
def MeshtoVTK(m,path):
    assert(isinstance(m,Mesh))
    points=[]
    polygons=[]
    for i in m.nodes:
        info = m.nodes[i]
        points.append(list(info))
    for i in m.elems.keys():
        n = m.elems[i][POS.TYPE]
        nlist = list(m.elems[i][POS.NODLIST])
        if(n==6):
            tmp=copy(nlist)
            nlist[2-1]=tmp[4-1]
            nlist[3-1]=tmp[2-1]
            nlist[4-1]=tmp[5-1]
            nlist[5-1]=tmp[3-1]       
        for j in range(len(nlist)):
            nlist[j]=nlist[j]-1
        polygons.append(nlist)
    s=PolyData(points=points,polygons=polygons)
    n=[]

#     for i in m.nodes:
        # info = m.nodes[i]
        # n.append(list(info[1:]))
    # pointdata=PointData(Normals(n,name='normal'))
    # vtk=VtkData(s,pointdata)
    vtk=VtkData(s)
    vtk.tofile(path,'ascii')

def ValuetoVTK(m,name,value):
    assert(isinstance(value,dict))
    points=[]
    polygons=[]
    for i in m.nodes:
        info = list(m.nodes[i])
        # replace z value with wave elevation
        info[2] = value[i] 
        points.append(info)

    for i in m.elems.keys():
        n = m.elems[i][POS.TYPE]
        nlist = list(m.elems[i][POS.NODLIST])
        if(n==6):
            tmp=copy(nlist)
            nlist[2-1]=tmp[4-1]
            nlist[3-1]=tmp[2-1]
            nlist[4-1]=tmp[5-1]
            nlist[5-1]=tmp[3-1]    
        for j in range(len(nlist)):
            nlist[j]=nlist[j]-1
        polygons.append(nlist)
    s=PolyData(points=points,polygons=polygons)
    vtk=VtkData(s)
    vtk.tofile(name,'ascii')


def waveVtk(m,path,fi):
    """
        waveVtk(m,path,fi)

        m: mesh
        path : output vtk name
        fi: input value file
    """
    assert(isinstance(m,Mesh))
    # assert(len(m.nodes)==len(value))
    with open(fi,"rb") as f:
        r = f.readlines()
        r1={}
        assert(len(m.nodes)==len(r))
        for i in range(len(r)):
            tmp=[float(j) for j in r[i].split()]
            #tmp[0] is node id
            r1[i+1]=tmp[1]

    ValuetoVTK(m,path,r1)
