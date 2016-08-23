import logging
from dutwav.__mesh_core import *
from dutwav.__mesh_core import _Mesh_core
from dutwav.draw import DrawMesh
from numpy import sqrt,round


class Mesh(_Mesh_core):
   def __init__(self,dp=4):
        _Mesh_core.__init__(self,dp)
        # self._num_fs_elem = 0
        # self._num_bd_elem = 0
        self.damp_info={}
        self.__rdrawObj=DrawMesh(self)
        # self._zonelist=[]
        # self._waterline=[]

        # self.xyz=None
        # self.dxyz=None
        # self.ncn=None

        self._waterline={}

   #
   def test_dp(self):
       print self.__cur_avail_nd_id

   # @func : update nrml vecotr for same tag group
   def _update_tag_nrml(self,taglist,vector):
      assert(len(vector)==3)
      nset=set()
      for i in self.elems:
          if self.elems[i][POS.TAG] in taglist:
              info=self.elems[i][POS.NRMLIST]
              for j in info:
                  nset.add(j)
              print nset    
      for i in nset:
          info = list(self.nrmls[i][0:1])+vector
          print info
          self.nrmls[i] = info 
      self._rebuild_rev_info()   

   # @func : reverse nrml direction for whole model
   def _reverse_nrml(self):
       for i in self.nrmls:
           info=list(self.nrmls[i])
           for j in [1,2,3]:
               info[j] = -info[j]
           self.nrmls[i]=tuple(info)

   # @func : give surface elems tag name 'surface'
   def _mark_surface_elems(self):
       self._mark_elems_at_z(0.0,'surface')

   # @func : mark element [pos]th position around value [v] with name [tag]
   def _mark_elems_at(self,pos,v,tag='marked'):
       s_node = set()
       s_elem = set()

       for i in self.nodes:
           if(abs(self.nodes[i][pos]-v)<1e-5):
               s_node.add(i)
       for i in self.elems:
           nodelist = set(self.elems[i][POS.NODLIST]) 
           if nodelist.issubset(s_node):
               s_elem.add(i)
       for i in s_elem:
           self.elems[i][0] = tag

   # @func : mark element with veritcal pos around [z] using name [tag]
   def _mark_elems_at_z(self,z,tag='marked'):
       s_node = set()
       s_elem = set()

       for i in self.nodes:
           if(abs(self.nodes[i][2]-z)<1e-5):
               s_node.add(i)
       for i in self.elems:
           nodelist = set(self.elems[i][POS.NODLIST]) 
           if nodelist.issubset(s_node):
               s_elem.add(i)
       for i in s_elem:
           self.elems[i][0] = tag

   # @func : list all nodes btw radius [rlow] and [rhigh]
   def _list_node_withr(self,rlow=0.,rhigh=1.):
       surface_node=set()
       for i in self.nodes:
           r = sqrt((self.nodes[i][1])**2+(self.nodes[i][0])**2)
           if(r<rhigh) and (r>rlow):
               surface_node.add(i)
       return(list(surface_node))   


   # @func : mark elems btw radius [rlow] and [rhigh]
   def _mark_elem_withr(self,rlow=0,rhigh=1,tag='marked'):
       surface_node = set(self._list_node_withr(rlow,rhigh))
       surface_elem = set()
       # for i in self.nodes:
           # r = sqrt((self.nodes[i][1])**2+(self.nodes[i][0])**2)
           # if(r<rhigh) and (r>rlow):
               # surface_node.add(i)
       for i in self.elems:
           nodelist = set(self.elems[i][POS.NODLIST]) 
           if nodelist.issubset(surface_node):
               surface_elem.add(i)

       for i in surface_elem:
           self.elems[i][0] = tag

   # @func: change groups in list [tag_list] to name [new_tag]
   def _update_tag(self,tag_list,new_tag):
       assert(isinstance(tag_list,list))
       assert(isinstance(new_tag,str))
       for i in self.elems:
           if self.elems[i][0] in tag_list:
               self.elems[i][0] = new_tag


   # @func: count elem number in each tag group
   def _count_elem(self):
       result={}
       for i in self.elems:
           key = self.elems[i][0]
           result[key] = result.get(key,0)+1
       return result 

   #===================================================
   #===================================================

   # @func: given damp function [f], define damp value for the surface nodes
   def _generate_damp_info(self,f=None):
       surface_node = set()
       for i in self.nodes:
          if(abs(self.nodes[i][2])<1e-5):
               surface_node.add(i)
       for i in surface_node:
          r = sqrt(self.nodes[i][0]**2+self.nodes[i][1]**2)
          if f:
              self.damp_info[i] = f(r)
          if not f:
              self.damp_info[i] = 0.0


   # @func: get set of waterline nodes
   def _get_waterline_node(self):
       print "elements in mesh must have well-defined tag\n\
               only 'body' elements will be looked into"
       nodes=set()
       for i in self.elems:
           info = self.elems[i]
           if info[0]=='body':
               nodelist = info[POS.NODLIST]
               for j in range(info[POS.TYPE]):
                   if self._is_node_fs(nodelist[j]):
                       nodes.add(nodelist[j])
       self._init_waterline(nodes)
       return(nodes)                

   def _get_btm_node(self,z=-1):
       print "elements in mesh must have well-defined tag\n\
               only 'body' elements will be looked into"
       nodes=set()
       for i in self.elems:
           info = self.elems[i]
           if info[0]=='body':
               nodelist = info[POS.NODLIST]
               for j in range(info[POS.TYPE]):
                   xyz=self.nodes[nodelist[j]]
                   if (abs(xyz[2]-z)<1e-5):
                       nodes.add(nodelist[j])
       return(nodes)                


   def _init_waterline(self,wset):
       for i in wset:
           key=round(list(self.nodes[i]),self.get_precision())
           self._waterline[tuple(key)] = None


   #===================================================
   #===================================================

   def validate_mesh(self):
      """
      first n node must be surface node
      first elems must be surface elems
      check coincident nodes
      check coincident nrmls
      check coincident elems
      """
      pass 

  #===================================================
   #===================================================

   def draw_model(self,points=[]):
       self.__rdrawObj.draw_model(points=points)

   def draw_lines(self,p=[],points=[]):
       self.__rdrawObj.draw_lines(p,points=points)

   def tecplt_value(self,path,value,soltime=1,kind =1):
       # if kind==2:
           # self.__rdrawObj.tecplt_value_poly(path,value)
       if kind==1:
           self.__rdrawObj.tecplt_value_quad(path,value,soltime)
       if kind==2:
           self.__rdrawObj.tecplt_value_poly_2(path,value)

   def tecplt_nrml(self,path,kind =1):
       if kind==2:
           # quad, node numbering, has nrml info
           self.__rdrawObj.tecplt_quad(path)
       if kind==1:
           #use nrml numbering, has nrml output
           self.__rdrawObj.tecplt_poly_3(path)
       if kind==3:
           #use node numbering, no nrml output
           self.__rdrawObj.tecplt_poly_2(path)
#        if kind==4:
           # #use node numbering, no nrml output
           # self.__rdrawObj.tecplt_poly_3(path)



   #===================================================
   #===================================================
   def extract_mesh(self,criteria):
       """
        new_mesh = extract_mesh([taglist])
        extract a new mesh from elements marked in taglist
       """
       assert(isinstance(criteria,list))
       # super class private cannot be accessed directly
       n = Mesh(self.get_precision())
       s_elem=set()
       s_node=set()
       s_nrml=set()
       # Gather all nrmls,nodes and elems from critia
       for e in self.elems:
           if self.elems[e][POS.TAG] in criteria:
               s_elem.add(e)
               s_node=s_node.union(set(self.elems[e][POS.NODLIST]))
               s_nrml=s_nrml.union(set(self.elems[e][POS.NRMLIST]))

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

       n._rebuild_rev_info()
       n._recreate_avail_info()

       return n


