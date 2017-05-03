from numpy import *
from dutwav.mesh import Mesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm


# input total number of elements on circle
n1 = 10
# input total number of elements along radius
n2 = 2

# input openning angle 
theta = 30
#input elements vertical
n3 = 2

# input depth 
# inner radius
# outer radius
depth =3
irad = 5
orad = 2



start_angle=0.+theta/2.
end_angle=360.-theta/2.

#points on circular
p1=deg2rad(linspace(start_angle,end_angle,2*n1+1))
# print rad2deg(p1)
#points on radical
p2=linspace(irad,orad,2*n2+1)
#points on vertical
p3=linspace(0,0-depth,2*n3+1)

m = Mesh()
dp = m.get_precision()

t1 = cos(p1)
t2 = sin(p1)

x=list()
y=list()
z=list()


#-------- Part 1--------------------
pms=(2*n1+1,2*n3+1)
ids = 1
#Generate the nodes
for j in range(pms[1]):
    for i in range(pms[0]):
        tmp=(round(irad*t1[i],dp),round(irad*t2[i],dp),round(p3[j],dp))
        m.nodes[ids]=tmp
        tmp=tmp/norm(tmp)
        m.nrmls[ids]=(ids,tmp[0],tmp[1],0)
        ids+=1

# Numbering the elements
ids = 1
for k in range(1,pms[1]-1,2):
    for j in range(1,pms[0]-1,2):
        i=j+(k-1)*pms[0]
        nodelist = (i,i+1,i+2,i+pms[0]+2,i+2*pms[0]+2,i+ \
                2*pms[0]+1,i+2*pms[0],i+pms[0])
        m.elems[ids] = ['mesh',8,nodelist,nodelist]
        ids+=1
m._recreate_avail_info()

#-------- Part 2--------------------
m1=Mesh()
ids=1
pms=(2*n1+1,2*n3+1)
#Generate the nodes
for j in range(pms[1]):
    for i in range(pms[0]):
        tmp=(round(orad*t1[i],dp),round(orad*t2[i],dp),round(p3[j],dp))
        m1.nodes[ids]=tmp
        tmp=tmp/norm(tmp)
        m1.nrmls[ids]=(ids,-tmp[0],-tmp[1],0)
        ids+=1
# Numbering the elements
ids = 1
for k in range(1,pms[1]-1,2):
    for j in range(1,pms[0]-1,2):
        i=j+(k-1)*pms[0]
        nodelist = (i,i+1,i+2,i+pms[0]+2,i+2*pms[0]+2,i+ \
                2*pms[0]+1,i+2*pms[0],i+pms[0])
        m1.elems[ids] = ['mesh',8,nodelist,nodelist]
        ids+=1
m.devour_mesh(m1)

#-------- Part 3--------------------
m1=Mesh()
ids=1
pms=(2*n1+1,2*n2+1)
#Generate the nodes
for j in range(pms[1]):
    for i in range(pms[0]):
        tmp=(round(p2[j]*t1[i],dp),round(p2[j]*t2[i],dp),round(p3[0],dp))
        m1.nodes[ids]=tmp
        m1.nrmls[ids]=(ids,0.,0.,1.)
        ids+=1
# Numbering the elements
ids = 1
for k in range(1,pms[1]-1,2):
    for j in range(1,pms[0]-1,2):
        i=j+(k-1)*pms[0]
        nodelist = (i,i+1,i+2,i+pms[0]+2,i+2*pms[0]+2,i+ \
                2*pms[0]+1,i+2*pms[0],i+pms[0])
        m1.elems[ids] = ['mesh',8,nodelist,nodelist]
        ids+=1

m.devour_mesh(m1)

#-------- Part 4--------------------
m1=Mesh()
ids=1
pms=(2*n1+1,2*n2+1)
#Generate the nodes
for j in range(pms[1]):
    for i in range(pms[0]):
        tmp=(round(p2[j]*t1[i],dp),round(p2[j]*t2[i],dp),round(p3[-1],dp))
        m1.nodes[ids]=tmp
        m1.nrmls[ids]=(ids,0.,0.,-1.)
        ids+=1
# Numbering the elements
ids = 1
for k in range(1,pms[1]-1,2):
    for j in range(1,pms[0]-1,2):
        i=j+(k-1)*pms[0]
        nodelist = (i,i+1,i+2,i+pms[0]+2,i+2*pms[0]+2,i+ \
                2*pms[0]+1,i+2*pms[0],i+pms[0])
        m1.elems[ids] = ['mesh',8,nodelist,nodelist]
        ids+=1

m.devour_mesh(m1)

#-------- Part 5--------------------
m1=Mesh()
ids=1
pms=(2*n2+1,2*n3+1)
#calc norm direction
alpha = float(theta)/2+90
# print "alpha=",alpha
alpha = deg2rad(alpha)
#Generate the nodes
for j in range(pms[1]):
    for i in range(pms[0]):
        tmp=(round(p2[i]*t1[0],dp),round(p2[i]*t2[0],dp),round(p3[j],dp))
        m1.nodes[ids]=tmp
        m1.nrmls[ids]=(ids,cos(alpha),sin(alpha),0.)
        ids+=1
# Numbering the elements
ids = 1
for k in range(1,pms[1]-1,2):
    for j in range(1,pms[0]-1,2):
        i=j+(k-1)*pms[0]
        nodelist = (i,i+1,i+2,i+pms[0]+2,i+2*pms[0]+2,i+ \
                2*pms[0]+1,i+2*pms[0],i+pms[0])
        m1.elems[ids] = ['mesh',8,nodelist,nodelist]
        ids+=1

m.devour_mesh(m1)

#----------------

#-------- Part 5--------------------
m1=Mesh()
ids=1
pms=(2*n2+1,2*n3+1)
#calc norm direction
alpha = -float(theta)/2-90
# print "alpha=",alpha
alpha = deg2rad(alpha)
#Generate the nodes
for j in range(pms[1]):
    for i in range(pms[0]):
        tmp=(round(p2[i]*t1[-1],dp),round(p2[i]*t2[-1],dp),round(p3[j],dp))
        m1.nodes[ids]=tmp
        m1.nrmls[ids]=(ids,cos(alpha),sin(alpha),0.)
        ids+=1
# Numbering the elements
ids = 1
for k in range(1,pms[1]-1,2):
    for j in range(1,pms[0]-1,2):
        i=j+(k-1)*pms[0]
        nodelist = (i,i+1,i+2,i+pms[0]+2,i+2*pms[0]+2,i+ \
                2*pms[0]+1,i+2*pms[0],i+pms[0])
        m1.elems[ids] = ['mesh',8,nodelist,nodelist]
        ids+=1

m.devour_mesh(m1)

#----------------u
#TODO remove unused nodes,and renumber model
# m1.nodes.pop(7)
print 'nodes=',len(m.nodes)
print 'nrmls=',len(m.nrmls)
m._renum_model()
m.tecplt_nrml("./ml.dat")

