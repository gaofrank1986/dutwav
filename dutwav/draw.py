import dutwav.mesh
import logging
from dutwav._mesh_core import POS
import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
        
class DrawMesh(object):
    def __init__(self,mesh):
        assert(isinstance(mesh,dutwav.mesh.Mesh))
        self.mesh = mesh
        print "only reference created,not copy"

    # This class is for drawing mesh func.



    def _draw_element(self,id,ax,scale,with_arrow,c='red'):
        this_elem = self.mesh.elems[id]
        loc = np.zeros((3,this_elem[POS.TYPE]+1),dtype='float')
        nrm = np.zeros((6,this_elem[POS.TYPE]),dtype='float')
        for i in range(this_elem[POS.TYPE]):
            loc[0:3,i] = self.mesh.nodes[this_elem[POS.NODLIST][i]][0:3]
            loc[:,-1] = loc[:,0] 
        ax.plot_wireframe(loc[0,:],loc[1,:],loc[2,:],rstride=5,cstride=5,color=c)
        if (with_arrow):
            for i in range(this_elem[POS.TYPE]):
                logging.debug(this_elem[POS.NRMLIST][i])
                nrm[0:3,i] = self.mesh.nrmls[this_elem[POS.NRMLIST][i]][1:4]
                nrm[3:6,i] = self.mesh.nodes[self.mesh.nrmls[this_elem[POS.NRMLIST][i]][0]][0:3]
                b = nrm[3:6,i]#base
                e = nrm[0:3,i]*scale+b#end
                a = Arrow3D([b[0],e[0]],[b[1],e[1]],[b[2],e[2]], mutation_scale=10, lw=1, arrowstyle="-|>", color="k")
                ax.add_artist(a)

    def _draw_element_damp(self,id,ax,scale,c='red'):
        this_elem = self.mesh.elems[id]
        loc = np.zeros((3,this_elem[POS.TYPE]+1),dtype='float')
        nrm = np.zeros((6,this_elem[POS.TYPE]),dtype='float')
        for i in range(this_elem[POS.TYPE]):
            loc[0:2,i] = self.mesh.nodes[this_elem[POS.NODLIST][i]][0:2]
            loc[2,i] = self.mesh.damp_info[this_elem[POS.NODLIST][i]]
            loc[:,-1] = loc[:,0] 
        ax.plot_wireframe(loc[0,:],loc[1,:],loc[2,:],rstride=5,cstride=5,color=c)
    
       
    def draw_model(self,scale=0.04,with_arrow=False,points=[],d = False):
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt
        from dutwav.util import set_aspect_equal_3d


        fig = plt.figure()
        # ax = fig.add_subplot(111,projection='3d')
        ax = fig.gca(projection='3d')
        # ax.set_aspect("equal")
        for i in self.mesh.elems:
            if not d:
                self._draw_element(i,ax,scale,with_arrow)
            else:
                self._draw_element_damp(i,ax,scale)
        # plt.show()
        if len(points)!=0:
            loc=np.zeros((3,len(points)))
            p=0
            for i in points:
                loc[0,p] = self.mesh.nodes[i][0]
                loc[1,p] = self.mesh.nodes[i][1]
                loc[2,p] = self.mesh.nodes[i][2]
                p+=1
            ax.scatter(loc[0,:],loc[1,:],loc[2,:],color="g",s=100)
        set_aspect_equal_3d(ax)
        plt.show()
        
        # return plt
    #==================================================================
    def export_tecplt_quad(self,path):
       with open(path,"wb") as f:
           num_pts = len(self.mesh.nrmls)
           num_elms = len(self.mesh.elems)
           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z","DX","DY","DZ"
                """)
           f.write('ZONE T="MESH" N=        {:d} ,E=         {:d} \
                   ,F=FEPOINT,ET=QUADRILATERAL\n'.format(num_pts,num_elms))
           for i in self.mesh.nrmls:
               info = self.mesh.nrmls[i]
               node = self.mesh.nodes[info[0]]
               f.write('   '.join('{0:<9.6f}'.format(j) for j in node[0:3]))
               f.write('   ')
               f.write('   '.join('{0:<9.6f}'.format(j) for j in info[1:4]))
               f.write("\n")
           for i in self.mesh.elems:
               info = self.mesh.elems[i][POS.NRMLIST]
               # nlist = [info[0],info[2],info[4],info[6]] 
               nlist = list(info)[::2]
               f.write('   '.join('{0:<8d}'.format(j) for j in nlist))
               f.write("\n")

    def export_tecplt_poly(self,path):
       with open(path,"wb") as f:
           num_pts = len(self.mesh.nrmls)
           num_elem = len(self.mesh.elems)
           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z","DX","DY","DZ"
                """)
           f.write('ZONE T="Mesh", ZONETYPE=FEPOLYGON, NODES= {:6d}, ELEMENTS= {:6d}, Faces= {:6d}, NumConnectedBoundaryFaces=0,TotalNumBoundaryConnections=0\n'.format(num_pts,num_elem,8*num_elem))
           # f.write('ZONE T="MESH" N=        {:d} ,E=         {:d} \
                   # ,F=FEPOINT,ET=QUADRILATERAL\n'.format(num_pts,num_elms))
           psl = []
           for i in range(8): 
               psl.append([])
           for i in self.mesh.nrmls:
               info = self.mesh.nrmls[i]
               node = self.mesh.nodes[info[0]]
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
           # str3=''
           # str2=''
           
           for i in self.mesh.elems.keys():
               n = self.mesh.elems[i][POS.TYPE]
               nlist = list(self.mesh.elems[i][POS.NRMLIST])
               nlist.append(nlist[0])
               print nlist
               print n
               for k in range(n):
                   # a = nlist[k]
                   # b = nlist[k+1]
                   # print a,b 
                   # f.write(' {0:<d} {0:<d}'.format(a,b))
                   f.write(' '.join('{0:<d}'.format(j) for j in nlist[k:k+2]))
                   f.write("\n")
                   psl[6].append(i)
                   psl[7].append(0)
                   # str3+=' {0:<d}'.format(int(i))
                   # str2+=' {0:<d}'.format(0)
               # if i%6 ==0 :
                   # str3+='\n'
               # if i%10 ==0 :
                   # str2+='\n'
           # str3+='\n'
           # str2+='\n'
           # f.write(str3)
           # f.write(str2)
           for i in [6,7]:
               max_len = len(psl[6])
               cha = max_len/500 + 1
               for k in range(cha):
                   f.write('   '.join('{0:<d}'.format(j) for j in psl[i][k*500:(k+1)*500]))
                   f.write('\n')


