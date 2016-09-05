import dutwav.mesh
from dutwav.__mesh_core import POS
# from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from numpy import mean,zeros

def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""

    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])
    
# class Arrow3D(FancyArrowPatch):
    # def __init__(self, xs, ys, zs, *args, **kwargs):
        # FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        # self._verts3d = xs, ys, zs

    # def draw(self, renderer):
        # xs3d, ys3d, zs3d = self._verts3d
        # xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        # self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        # FancyArrowPatch.draw(self, renderer)
        
class DrawMesh(object):
    def __init__(self,mesh):
        assert(isinstance(mesh,dutwav.mesh.Mesh))
        self.__rMeshObj = mesh

    # This class is for drawing mesh func.


    def _draw_nrml(self,id,ax,scale):
        nrm[0:3,i] = self.__rMeshObj.nrmls[id][1:4]
        nrm[3:6,i] = self.__rMeshObj.nodes[self.__rMeshObj.nrmls[id][0]][0:3]
        b = nrm[3:6,i]#base
        e = nrm[0:3,i]*scale+b#end
        # a = Arrow3D([b[0],e[0]],[b[1],e[1]],[b[2],e[2]], mutation_scale=10, lw=1, arrowstyle="-|>", color="k")
        # ax.add_artist(a)

        
    def _draw_element(self,id,ax,c='red'):
        this_elem = self.__rMeshObj.elems[id]
        loc = zeros((3,this_elem[POS.TYPE]+1),dtype='float')
        nrm = zeros((6,this_elem[POS.TYPE]),dtype='float')
        for i in range(this_elem[POS.TYPE]):
            loc[0:3,i] = self.__rMeshObj.nodes[this_elem[POS.NODLIST][i]][0:3]
            loc[:,-1] = loc[:,0] 
        ax.plot_wireframe(loc[0,:],loc[1,:],loc[2,:],rstride=5,cstride=5,color=c)
           
    def draw_lines(self,pi=[],s=10,points=[]):

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for i in range(len(pi)-1):
            n1=pi[i]
            n2=pi[i+1]
            print n1,n2
            self._draw_line(n1,n2,ax,s=s)
        if len(points)!=0:
            loc=zeros((3,len(points)))
            p=0
            for i in points:
                loc[0,p] = self.__rMeshObj.nodes[i][0]
                loc[1,p] = self.__rMeshObj.nodes[i][1]
                loc[2,p] = self.__rMeshObj.nodes[i][2]
                p+=1
            ax.scatter(loc[0,:],loc[1,:],loc[2,:],color="g",s=s*5)

        set_aspect_equal_3d(ax)
        plt.show()

    def _draw_line(self,n1,n2,ax,c='red',s=10):
        info1 = list(self.__rMeshObj.nodes[n1])
        info2 = list(self.__rMeshObj.nodes[n2])
        loc = zeros((3,2),dtype='float')
        loc[0:3,0] = info1[0:3]
        loc[0:3,1] = info2[0:3]
        ax.plot_wireframe(loc[0,:],loc[1,:],loc[2,:],rstride=5,cstride=5,color=c)          
        ax.scatter(loc[0,:],loc[1,:],loc[2,:],color="g",s=s)

    # def _draw_element_damp(self,id,ax,scale,c='red'):
        # this_elem = self.__rMeshObj.elems[id]
        # loc = zeros((3,this_elem[POS.TYPE]+1),dtype='float')
        # nrm = zeros((6,this_elem[POS.TYPE]),dtype='float')
        # for i in range(this_elem[POS.TYPE]):
            # loc[0:2,i] = self.__rMeshObj.nodes[this_elem[POS.NODLIST][i]][0:2]
            # loc[2,i] = self.__rMeshObj.damp_info[this_elem[POS.NODLIST][i]]
            # loc[:,-1] = loc[:,0] 
        # ax.plot_wireframe(loc[0,:],loc[1,:],loc[2,:],rstride=5,cstride=5,color=c)
    
       
    def draw_model(self,scale=0.04,with_arrow=False,points=[],d = False):


        fig = plt.figure()
        # ax = fig.add_subplot(111,projection='3d')
        ax = fig.gca(projection='3d')
        # ax.set_aspect("equal")
        for i in self.__rMeshObj.elems:
            if not d:
                self._draw_element(i,ax)
            else:
                self._draw_element_damp(i,ax)
        if with_arrow:
            for i in self.nrmls:
                _draw_nrml(i,ax,scale)
        if len(points)!=0:
            loc=zeros((3,len(points)))
            p=0
            for i in points:
                loc[0,p] = self.__rMeshObj.nodes[i][0]
                loc[1,p] = self.__rMeshObj.nodes[i][1]
                loc[2,p] = self.__rMeshObj.nodes[i][2]
                p+=1
            ax.scatter(loc[0,:],loc[1,:],loc[2,:],color="g",s=100)
        set_aspect_equal_3d(ax)
        plt.show()
        
    #==================================================================

    # @func : export tecplot mesh using quad element
    def tecplt_quad(self,path):
       with open(path,"wb") as f:
           num_pts = len(self.__rMeshObj.nrmls)
           num_elms = len(self.__rMeshObj.elems)
           f.write("""TITLE = "3D Mesh Grid Data for Element Boundary"\n\
                   VARIABLES = "X", "Y", "Z","DX","DY","DZ" \n""")
           f.write('ZONE T="MESH" N=        {:d} ,E=         {:d} \
                   ,F=FEPOINT,ET=QUADRILATERAL\n'.format(num_pts,num_elms))
       
           for i in self.__rMeshObj.nrmls:
               info = self.__rMeshObj.nrmls[i]
               node = self.__rMeshObj.nodes[info[0]]
               f.write('   '.join('{0:<9.6f}'.format(j) for j in node[0:3]))
               f.write('   ')
               f.write('   '.join('{0:<9.6f}'.format(j) for j in info[1:4]))
               f.write("\n")
           for i in self.__rMeshObj.elems:
               n = self.__rMeshObj.elems[i][POS.TYPE]
               info = list(self.__rMeshObj.elems[i][POS.NRMLIST])
               # nlist = [info[0],info[2],info[4],info[6]] 
               nlist = info[::2]
               if(n==6):
                   nlist=[info[0],info[1],info[2],info[0]]
               f.write('   '.join('{0:<8d}'.format(j) for j in nlist))
               f.write("\n")



    # @func : export surface mesh using quad, also adding time info
    def tecplt_value_quad(self,path,value,soltime=1):

       assert(isinstance(value,list))
       for i in value:
           assert(len(self.__rMeshObj.nodes)==len(i))

       with open(path,"wb") as f:
           num_pts = len(self.__rMeshObj.nodes)
           num_elms = len(self.__rMeshObj.elems)
           str1='TITLE = "3D Mesh Grid Data for Element Boundary"\n'+ \
                  'VARIABLES = "X", "Y","Z"'
           for i in range(len(value)):
               str1+=', "v'+str(i+1)+'"'
           str1+='\n'
           f.write(str1)
           f.write('ZONE T="MESH" N=    {:d} ,E=     {:d}'.format(num_pts,num_elms) \
                   +',F=FEPOINT,ET=QUADRILATERAL\n')
           f.write('SolutionTime={:d}\n'.format(soltime))
           f.write('StrandID={:d}\n'.format(1))

           for i in self.__rMeshObj.nodes:
               node = self.__rMeshObj.nodes[i]
               # output x,y,z
               f.write('   '.join('{0:<9.6f}'.format(j) for j in node[0:3]))
               # output value info
               for j in range(len(value)):
                   f.write('    ')
                   f.write('{0:<9.6f}'.format(value[j][i]))
               f.write("\n")
            # output node connection list
           for i in self.__rMeshObj.elems:
               info = list(self.__rMeshObj.elems[i][POS.NODLIST])
               typ = self.__rMeshObj.elems[i][POS.TYPE]
               nlist = info[::2]
               if (typ==6):
                   nlist=[info[0],info[1],info[2],info[0]]
               f.write('   '.join('{0:<8d}'.format(j) for j in nlist))
               f.write("\n")





    # @func: use polygon elem,have nrml info,no time info
    # plot tecplot, use node numbering
    def tecplt_poly_2(self,path):
       print "Nodes Numbering Used"
       print "No nrml info,No time output"
       MAX_LINE_HOLD=500
       with open(path,"wb") as f:
           num_pts = len(self.__rMeshObj.nodes)
           num_elem = len(self.__rMeshObj.elems)
           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z"\n
                """)
           f.write('ZONE T="Mesh", ZONETYPE=FEPOLYGON, NODES= {:6d}, ELEMENTS= {:6d}, Faces= {:6d}, NumConnectedBoundaryFaces=0,TotalNumBoundaryConnections=0\n'.format(num_pts,num_elem,8*num_elem))
           
           psl = []
           for i in range(8): 
               psl.append([])

           for i in self.__rMeshObj.nodes:
               node = self.__rMeshObj.nodes[i]
               for j in [0,1,2]:
                   psl[j].append(node[j])

           for i in range(3):#3 is number of psl to output
               max_len = len(psl[0])
               cha = max_len/MAX_LINE_HOLD + 1
               for k in range(cha):
                   f.write('   '.join('{0:<7.4f}'.format(j) for j in psl[i][k*MAX_LINE_HOLD:(k+1)*MAX_LINE_HOLD]))
                   f.write('\n')
           
           for i in self.__rMeshObj.elems.keys():
               n = self.__rMeshObj.elems[i][POS.TYPE]
               nlist = list(self.__rMeshObj.elems[i][POS.NODLIST])
               nlist.append(nlist[0])
               for k in range(n):
                   f.write(' '.join('{0:<d}'.format(j) for j in nlist[k:k+2]))
                   f.write("\n")
                   psl[6].append(i)
                   psl[7].append(0)

           for i in [6,7]:
               max_len = len(psl[6])
               cha = max_len/MAX_LINE_HOLD + 1
               for k in range(cha):
                   f.write('   '.join('{0:<d}'.format(j) for j in psl[i][k*MAX_LINE_HOLD:(k+1)*MAX_LINE_HOLD]))
                   f.write('\n')


    # @func : export tecplot mesh using polygon element
    # use nrml for numbering
    def tecplt_poly_1(self,path):
       print "poly_1"
       print "nrml Numbering Used"
       print "No time output"
       with open(path,"wb") as f:
           num_pts = len(self.__rMeshObj.nrmls)
           num_elem = len(self.__rMeshObj.elems)
           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z","DX","DY","DZ"\n
                """)
           f.write('ZONE T="Mesh", ZONETYPE=FEPOLYGON, NODES= {:6d}, ELEMENTS= {:6d}, Faces= {:6d}, NumConnectedBoundaryFaces=0,TotalNumBoundaryConnections=0\n'.format(num_pts,num_elem,8*num_elem))

           psl = []
           for i in range(8): 
               psl.append([])
           #use nrml id for numbering
           for i in self.__rMeshObj.nrmls:
               info = self.__rMeshObj.nrmls[i]
               node = self.__rMeshObj.nodes[info[0]]
               for j in [0,1,2]:
                   psl[j].append(node[j])
                   psl[j+3].append(info[j+1]) 

           # print psl[0] 
           for i in range(6):
               max_len = len(psl[0])
               cha = max_len/500 + 1
               for k in range(cha):
                   f.write('   '.join('{0:<7.4f}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   f.write('\n')

           
           for i in self.__rMeshObj.elems.keys():
               n = self.__rMeshObj.elems[i][POS.TYPE]
               nlist = list(self.__rMeshObj.elems[i][POS.NRMLIST])
               #add first node
               nlist.append(nlist[0])
               for k in range(n):
                   f.write(' '.join('{0:<d}'.format(j) for j in nlist[k:k+2]))
                   f.write("\n")
                   psl[6].append(i)
                   psl[7].append(0)

           for i in [6,7]:
               max_len = len(psl[6])
               cha = max_len/500 + 1
               for k in range(cha):
                   f.write('   '.join('{0:<d}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   f.write('\n')

    # @func : export surface mesh using polygon,no time info,no nrml info
    def tecplt_value_poly(self,path,value):
       print "Nodes Numbering Used"
       print "No Nrml Ouput,No time output"
       with open(path,"wb") as f:
           num_pts = len(self.__rMeshObj.nodes)
           num_elem = len(self.__rMeshObj.elems)
           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z","v1"\n
                """)
           f.write('ZONE T="Mesh", ZONETYPE=FEPOLYGON, NODES= {:6d}, ELEMENTS= {:6d}, Faces= {:6d}, NumConnectedBoundaryFaces=0,TotalNumBoundaryConnections=0\n'.format(num_pts,num_elem,8*num_elem))

           psl = []
           for i in range(8): 
               psl.append([])

           for i in self.__rMeshObj.nodes:
               node = self.__rMeshObj.nodes[i]
               for j in [0,1,2]:
                   psl[j].append(node[j])
               psl[3].append(value[0][i])

           # psl[1..3]
           for i in range(4):
               max_len = len(psl[0])
               cha = max_len/500 + 1
               for k in range(cha):
                   f.write('   '.join('{0:<7.4f}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   f.write('\n')
           
           for i in self.__rMeshObj.elems.keys():
               n = self.__rMeshObj.elems[i][POS.TYPE]
               nlist = list(self.__rMeshObj.elems[i][POS.NODLIST])
               nlist.append(nlist[0])
               for k in range(n):
                   f.write(' '.join('{0:<d}'.format(j) for j in nlist[k:k+2]))
                   f.write("\n")
                   psl[6].append(i)
                   psl[7].append(0)
           for i in [6,7]:
               max_len = len(psl[6])
               cha = max_len/500 + 1
               for k in range(cha):
                   f.write('   '.join('{0:<d}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   f.write('\n')




 # @func : export tecplot mesh using polygon element
    # use nrml for numbering
    def tecplt_poly_3(self,path):
       print "poly_3"
       print "nrml Numbering Used"
       print "No time output"
       from copy import copy
       with open(path,"wb") as f:
           num_pts = len(self.__rMeshObj.nrmls)
           num_elem = len(self.__rMeshObj.elems)


           bodystr=''
           psl = []
           for i in range(8): 
               psl.append([])
           #use nrml id for numbering
           for i in self.__rMeshObj.nrmls:
               info = self.__rMeshObj.nrmls[i]
               node = self.__rMeshObj.nodes[info[0]]
               for j in [0,1,2]:
                   psl[j].append(node[j])
                   psl[j+3].append(info[j+1]) 

           # print psl[0] 
           for i in range(6):
               max_len = len(psl[0])
               cha = max_len/500 + 1
               for k in range(cha):
                   bodystr += ('   '.join('{0:<7.4f}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   bodystr +=('\n')

           
           for i in self.__rMeshObj.elems.keys():
               n = self.__rMeshObj.elems[i][POS.TYPE]
               nlist = list(self.__rMeshObj.elems[i][POS.NRMLIST])
               if(n==6):
                   tmp=copy(nlist)
                   nlist[2-1]=tmp[4-1]
                   nlist[3-1]=tmp[2-1]
                   nlist[4-1]=tmp[5-1]
                   nlist[5-1]=tmp[3-1]

                    
               #add first node
               nlist.append(nlist[0])
               for k in range(n):
                   bodystr+=(' '.join('{0:<d}'.format(j) for j in nlist[k:k+2]))
                   bodystr+=("\n")
                   psl[6].append(i)
                   psl[7].append(0)

           for i in [6,7]:
               max_len = len(psl[6])
               cha = max_len/500 + 1
               for k in range(cha):
                   bodystr+=('   '.join('{0:<d}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   bodystr+=('\n')

           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z","DX","DY","DZ"\n
                """)
           f.write('ZONE T="Mesh", ZONETYPE=FEPOLYGON, NODES= {:6d}, ELEMENTS= {:6d}, Faces= {:6d}, NumConnectedBoundaryFaces=0,TotalNumBoundaryConnections=0\n'.format(num_pts,num_elem,len(psl[6])))
           f.write(bodystr)


               # @func : export surface mesh using polygon,no time info,no nrml info
    def tecplt_value_poly_2(self,path,value):
       print "Nodes Numbering Used"
       print "No Nrml Ouput,No time output"
       with open(path,"wb") as f:
           num_pts = len(self.__rMeshObj.nodes)
           num_elem = len(self.__rMeshObj.elems)
           ss=''

           psl = []
           for i in range(8): 
               psl.append([])

           for i in self.__rMeshObj.nodes:
               node = self.__rMeshObj.nodes[i]
               for j in [0,1,2]:
                   psl[j].append(node[j])
               psl[3].append(value[0][i])

           # psl[1..3]
           for i in range(4):
               max_len = len(psl[0])
               cha = max_len/500 + 1
               for k in range(cha):
                   ss+=('   '.join('{0:<7.4f}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   ss+=('\n')
           
           for i in self.__rMeshObj.elems.keys():
               n = self.__rMeshObj.elems[i][POS.TYPE]
               nlist = list(self.__rMeshObj.elems[i][POS.NODLIST])
               if(n==6):
                   tmp=copy(nlist)
                   nlist[2-1]=tmp[4-1]
                   nlist[3-1]=tmp[2-1]
                   nlist[4-1]=tmp[5-1]
                   nlist[5-1]=tmp[3-1]
               nlist.append(nlist[0])
               for k in range(n):
                   ss+=(' '.join('{0:<d}'.format(j) for j in nlist[k:k+2]))
                   ss+=("\n")
                   psl[6].append(i)
                   psl[7].append(0)
           for i in [6,7]:
               max_len = len(psl[6])
               cha = max_len/500 + 1
               for k in range(cha):
                   ss+=('   '.join('{0:<d}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   ss+=('\n')

           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z","v1"\n
                """)
           f.write('ZONE T="Mesh", ZONETYPE=FEPOLYGON, NODES= {:6d}, ELEMENTS= {:6d}, Faces= {:6d}, NumConnectedBoundaryFaces=0,TotalNumBoundaryConnections=0\n'.format(num_pts,num_elem,len(psl[6])))

           f.write(ss)
