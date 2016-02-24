import logging
import numpy as np
from dutwav._mesh_core import *
from dutwav._mesh_core import _Mesh_core


class Mesh(_Mesh_core):
   def __init__(self,dp=4):
        _Mesh_core.__init__(self,dp)
        # self._num_fs_elem = 0
        # self._num_bd_elem = 0
        self.damp_info={}
        # self._zonelist=[]
        # self._waterline=[]
        
        # self.xyz=None
        # self.dxyz=None
        # self.ncn=None

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
           key=np.round(list(self.nodes[i]),self.__dp)
           self._waterline[tuple(key)] = None

                    


   def test_mesh(self):
      """
      first n node must be surface node
      first elems must be surface elems
      """
      pass 



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
       n.__cur_avail_el_id = len(s_elem)+1
       n.__cur_avail_nd_id = len(s_node)+1
       n.__cur_avail_nm_id = len(s_nrml)+1
       # Rebuild rev_node,rev_rnml info
       # for i in n.nodes:
          # xyz = n.nodes[i]
          # key = (round(xyz[0],n.__dp),round(xyz[1],n.__dp),round(xyz[2],n.__dp))
          # n.rev_nd[key] = i
       # for i in n.nrmls:
          # info = n.nrmls[i]
          # key = (info[0],round(info[1],n.__dp),round(info[2],n.__dp),round(info[3],n.__dp))
          # n.rev_nm[key] = i
       n._rebuild_rev_info()
       
       return n

 
