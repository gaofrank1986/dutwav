import logging
import numpy as np
from dutwav._mesh_core import *
from dutwav._mesh_core import _Mesh_core

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

class Mesh(_Mesh_core):
   def __init__(self,dp=4):
        _Mesh_core.__init__(self,dp)
        self._num_fs_elem = 0
        self._num_bd_elem = 0
        self.damp_info={}
        # self._zonelist=[]
        # self._waterline=[]
        
        self.xyz=None
        self.dxyz=None
        self.ncn=None

        self._waterline={}
  
   def _mark_surface_elems(self):
       surface_node = set()
       surface_elem = set()
       
       for i in self.nodes:
           if(abs(self.nodes[i][2])<1e-5):
               surface_node.add(i)

       for i in self.elems:
           nodelist = set(self.elems[i][POS.NODLIST]) 
           if nodelist.issubset(surface_node):
               surface_elem.add(i)
       
       for i in surface_elem:
           self.elems[i][0] = "surface"

   def _mark_elem_withr(self,rlow=0,rhigh=1,tag='marked'):
       surface_node = set()
       surface_elem = set()
       
       for i in self.nodes:
           r = np.sqrt((self.nodes[i][1])**2+(self.nodes[i][0])**2)
           if(r<rhigh) and (r>rlow):
               surface_node.add(i)

       for i in self.elems:
           nodelist = set(self.elems[i][POS.NODLIST]) 
           if nodelist.issubset(surface_node):
               surface_elem.add(i)
       
       for i in surface_elem:
           self.elems[i][0] = tag

   def _update_tag(self,tag_list,new_tag):
      for i in self.elems:
         if self.elems[i][0] in tag_list:
            self.elems[i][0] = new_tag
        
   def _count_elem(self,tag=[]):
       result={}
       for i in self.elems:
           key = self.elems[i][0]
           result[key] = result.get(key,0)+1
       return result 
           
   def _generate_damp_info(self,f=None):
       surface_node = set()
       for i in self.nodes:
          if(abs(self.nodes[i][2])<1e-5):
             surface_node.add(i)
       for i in surface_node:
          r = np.sqrt(self.nodes[i][0]**2+self.nodes[i][1]**2)
          if f:
            self.damp_info[i] = f(r)
          if not f:
            self.damp_info[i] = 0.0



   def _get_waterline_node(self):
       print "elements in mesh must have well-defined tag\n\
       only 'body' elements will be looked into"
       nodes=set()
       for i in self.elems:
           info = self.elems[i]
           if info[0]=='body':
               nodelist = info[POS.NODLIST]
               for j in range(info[POS.ELEMTYPE]):
                   if self._is_node_fs(nodelist[j]):
                       nodes.add(nodelist[j])
       self._init_waterline(nodes)
       return(nodes)                

   def _init_waterline(self,wset):
       for i in wset:
           key=np.round(list(self.nodes[i]),self._dp)
           self._waterline[tuple(key)] = None

                    


   def _draw_element(self,id,ax,scale,with_arrow,c='red'):
      this_elem = self.elems[id]
      loc = np.zeros((3,this_elem[POS.TYPE]+1),dtype='float')
      nrm = np.zeros((6,this_elem[POS.TYPE]),dtype='float')
      for i in range(this_elem[POS.TYPE]):
         loc[0:3,i] = self.nodes[this_elem[POS.NODLIST][i]][0:3]
         loc[:,-1] = loc[:,0] 
      ax.plot_wireframe(loc[0,:],loc[1,:],loc[2,:],rstride=5,cstride=5,color=c)
      if (with_arrow):
         for i in range(this_elem[1]):
            logging.debug(this_elem[POS.NRMLIST][i])
            nrm[0:3,i] = self.nrmls[this_elem[POS.NRMLIST][i]][1:4]
            nrm[3:6,i] = self.nodes[self.nrmls[this_elem[POS.NRMLIST][i]][0]][0:3]
            b = nrm[3:6,i]#base
            e = nrm[0:3,i]*scale+b#end
            a = Arrow3D([b[0],e[0]],[b[1],e[1]],[b[2],e[2]], mutation_scale=10, lw=1, arrowstyle="-|>", color="k")
            ax.add_artist(a)

   def _draw_element_damp(self,id,ax,scale,c='red'):
      this_elem = self.elems[id]
      loc = np.zeros((3,this_elem[POS.TYPE]+1),dtype='float')
      nrm = np.zeros((6,this_elem[POS.TYPE]),dtype='float')
      for i in range(this_elem[POS.TYPE]):
         loc[0:2,i] = self.nodes[this_elem[POS.NODLIST][i]][0:2]
         loc[2,i] = self.damp_info[this_elem[POS.NODLIST][i]]
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
      for i in self.elems:
          if not d:
             self._draw_element(i,ax,scale,with_arrow)
          else:
             self._draw_element_damp(i,ax,scale)
      # plt.show()
      if len(points)!=0:
          loc=np.zeros((3,len(points)))
          p=0
          for i in points:
              loc[0,p] = self.nodes[i][0]
              loc[1,p] = self.nodes[i][1]
              loc[2,p] = self.nodes[i][2]
              p+=1
          ax.scatter(loc[0,:],loc[1,:],loc[2,:],color="g",s=100)
      set_aspect_equal_3d(ax)
      plt.show()
      # return plt


   def export_tecplt_quad(self,path):
       with open(path,"wb") as f:
           num_pts = len(self.nrmls)
           num_elms = len(self.elems)
           f.write("""
           TITLE = "3D Mesh Grid Data for Element Boundary"
            VARIABLES = "X", "Y", "Z","DX","DY","DZ"
                """)
           f.write('ZONE T="MESH" N=        {:d} ,E=         {:d} \
                   ,F=FEPOINT,ET=QUADRILATERAL\n'.format(num_pts,num_elms))
           for i in self.nrmls:
               info = self.nrmls[i]
               node = self.nodes[info[0]]
               f.write('   '.join('{0:<9.6f}'.format(j) for j in node[0:3]))
               f.write('   ')
               f.write('   '.join('{0:<9.6f}'.format(j) for j in info[1:4]))
               f.write("\n")
           for i in self.elems:
               info = self.elems[i][POS.NRMLIST]
               # nlist = [info[0],info[2],info[4],info[6]] 
               nlist = list(info)[::2]
               f.write('   '.join('{0:<8d}'.format(j) for j in nlist))
               f.write("\n")

   def export_tecplt_poly(self,path):
       with open(path,"wb") as f:
           num_pts = len(self.nrmls)
           num_elem = len(self.elems)
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
           for i in self.nrmls:
               info = self.nrmls[i]
               node = self.nodes[info[0]]
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
           
           for i in self.elems.keys():
               n = self.elems[i][POS.TYPE]
               nlist = list(self.elems[i][POS.NRMLIST])
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



   def output(self,path,kind):
      if kind=='b':
         self._output_as_bd(path)
      elif kind == 'f':
         self._output_as_fs(path)
      else:
         print "Please enter \'b\' for body format or \'f\' for free surface format"

               

   def shift_model(self,vector):
       for i in self.nodes:
           xyz = np.array(self.nodes[i],dtype='float64')
           xyz = np.array(vector,dtype='float64')+xyz
           self.nodes[i]=tuple(xyz)
   
   def _redo_bd_nrml(self):
      import scipy.linalg as la
      flag = True
      for i in self.elems:
         if self.elems[i][0] != 'body':
            flag = False
      if not flag:
         print "Please make sure all elements is defined as 'body'"
         return
      self.rev_nm={}
      for i in self.nrmls:
         info = self.nrmls[i]
         xyz = np.array(self.nodes[info[0]],dtype='float64')
         xyz[2] = 0.
         nrml = xyz/la.norm(xyz)#normalize vector 
         self.nrmls[i] = (info[0],nrml[0],nrml[1],nrml[2])
         nrml = np.round(nrml,self._dp)
         key = (info[0],nrml[0],nrml[1],nrml[2])
         self.rev_nm[key] = i
      return
         
   def _update_all_nrml(self,vector):
      print "vecot will be auto normalized"
      import scipy.linalg as la
      self.rev_nm={}
      for i in self.nrmls:
         info = self.nrmls[i]
         xyz = np.array(vector,dtype='float64')
         nrml = xyz/la.norm(xyz)#normalize vector 
         self.nrmls[i] = (info[0],nrml[0],nrml[1],nrml[2])
         nrml = np.round(nrml,self._dp)
         key = (info[0],nrml[0],nrml[1],nrml[2])
         self.rev_nm[key] = i

          


   def extract_mesh(self,criteria):
       """
        new_mesh = extract_mesh([taglist])
        extract a new mesh from elements marked in taglist
       """
            
       import dutwav.mesh as ms
       n = ms.Mesh()
       s_elem=set()
       s_node=set()
       s_nrml=set()
       # Gather all nrmls,nodes and elems from critia
       for e in self.elems:
           if self.elems[e][POS.TAG] in criteria:
               s_elem.add(e)
               s_node=s_node.union(set(self.elems[e][POS.NODLIST]))
               s_nrml=s_nrml.union(set(self.elems[e][POS.NRMLIST]))
       # 
       s_node = sorted(list(s_node))
       renum_node = {}
       for i in range(len(s_node)):
          n.nodes[i+1] = self.nodes[s_node[i]]
          renum_node[s_node[i]]=i+1
       # ============================   
       s_nrml = sorted(list(s_nrml))
       renum_nrml={}
       for i in range(len(s_nrml)):
          info = list(self.nrmls[s_nrml[i]])
          info[0] = renum_node[info[0]]
          n.nrmls[i+1] = info
          renum_nrml[s_nrml[i]]=i+1
       #
       s_elem = sorted(list(s_elem))
       for i in range(len(s_elem)):
           info = self.elems[s_elem[i]]
           nodelist = list(info[2])
           nrmlist = list(info[3])
           for j in range(info[1]):
              nodelist[j] = renum_node[nodelist[j]]
              nrmlist[j] = renum_nrml[nrmlist[j]]
           n.elems[i+1] = ['extract',info[1],tuple(nodelist),tuple(nrmlist)]
       n._cur_avail_el_id = len(s_elem)+1
       n._cur_avail_nd_id = len(s_node)+1
       n._cur_avail_nm_id = len(s_nrml)+1
       # Rebuild rev_node,rev_rnml info
       # for i in n.nodes:
          # xyz = n.nodes[i]
          # key = (round(xyz[0],n._dp),round(xyz[1],n._dp),round(xyz[2],n._dp))
          # n.rev_nd[key] = i
       # for i in n.nrmls:
          # info = n.nrmls[i]
          # key = (info[0],round(info[1],n._dp),round(info[2],n._dp),round(info[3],n._dp))
          # n.rev_nm[key] = i
       n._rebuild_rev_info()
       
       return n

   def mirror_mesh(self,kind=2,base=[0,0,0]):
       """mirror_mesh(kind=2,base=[0,0,0])
         kind=0,along yz plane
         kind=1,along xz plane
         kind=2,along xz plane
       """
       import copy
       n = copy.deepcopy(self)
       for i in n.nodes:
           info = list(n.nodes[i])
           info[kind] = 2*base[kind]-info[kind] 
           n.nodes[i]=tuple(info)
       for i in n.nrmls:
           info = list(n.nrmls[i])
           info[kind+1] =-info[kind+1] 
           n.nrmls[i]=tuple(info)
       n._rebuild_rev_info()
       return n

