import logging
from enum import Enum
from copy import deepcopy
from scipy.linalg import norm
import numpy as np
from numpy import cos,sin,deg2rad,array,sqrt
from numpy import round as npround


class POS(Enum):
    TAG=0
    TYPE=1
    NODLIST=2
    NRMLIST=3

class FILETYPE(Enum):
    BODY=0
    SURFACE_WO_DAMP=1
    SURFACE_W_DAMP=2
    SURFACE_AUTO_DAMP=3
    EXTERNAL = 9

# class Node(object):
    # # __slots__={'pos','id'}
    # def __init__(self,id,pos):
        # self.pos= pos
        # self.id = id
        # self.connected_elem=[]
        # self.connected_edge=[]
        # self.associated_nrml=[]
        # self.coor = None
    
    # def print_node(self):
        # pass

# class Normal(object):
    # # __slots__={'pos','id'}
    # def __init__(self,id,pos):
        # self.pos= pos
        # self.id = id
        # self.associated_node=[]
    
    # def print_nrml(self):
        # pass

# class Coordinate(object):
    # def __init__(self,id,ctype,offset):
        # self.id = id
        # self.type = ctype
        # self.offset = offset

    # def print_coor():
        # pass

# ====================================
# class Element(object):
    # # __slots__={'nl','nrl','etype','kind','id'}
   # def __init__(self,ename):
      # self.tag = tag
      # self.etype = -1#6 for triangle,8 for rectangular,4 for 4 node rect,3 for 3 node triangle
      # self.nl=[]
      # self.nrl=[]

   # def print_elem(self,flag=0): # flag =0  for nodelist,flag=1 for nrml list
        # pass


#=============================================
'''
    define the most basic data strucutre for mesh object
    Mesh---------
    nsys   : symmetric info /currently not used
    nodes  : [dict] contains tuple (x,y,z)
    nodes  : [dict] contains tuple (base,x,y,z)
    elems  : [dict] tuple('type',node_count,[nodelist],[nrmlist])
    rev_nd : [dict] given key (x,y,z) rounded for __dp, return node id
    rev_nm : refer to above
    edges  : [dict] tuple (node id 1,node id 2) 1 is smaller than 2,ordered 
'''

class _Mesh_core(object):
   def __init__(self,dp=4):
        self.nsys = 0

        self.nodes = {}
        self.elems = {}
        self.rev_nd = {}
        self.rev_nm = {}
        self.nrmls ={}
        # self.coors = {}
        self.edges={}

        self.__cur_avail_nd_id=1
        self.__cur_avail_el_id=1
        self.__cur_avail_nm_id=1
        self.__dp=dp

   def get_precision(self):
       '''
          __dp is used for rounding up, used in rev_nm,rev_nd
       '''
       return self.__dp

   def _is_exist_node(self,pos):
       '''
       '''
       assert(len(pos)==3)
       key = (round(pos[0],self.__dp),round(pos[1],self.__dp),round(pos[2],self.__dp))
       return(key in self.rev_nd)

   def _is_exist_nrml(self,node_id,pos):
       assert(len(pos)==3)
       key = (node_id,round(pos[0],self.__dp),round(pos[1],self.__dp),round(pos[2],self.__dp))
       return(key in self.rev_nd)

   #see if node is on free surface
   def _is_node_fs(self,n,z=0):
       xyz=self.nodes[n]
       if abs(xyz[2]-z)<1e-7:
           return True
       return False
   
        
   

   def _redo_bd_nrml(self):
      """
            Create nrml for cylinder body (without bottom) ,pointing outward
            -----------------------
            Note: only check those element with 'body' tag
            Note: all elements should be taged 'body'
      """
      print ("Create nrml for cylinder body (without bottom) ,pointing outward")
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
         # xyz = array(self.nodes[info[0]],dtype='float64')
         xyz=np.array(self.nodes[info[0]])
         xyz[2] = 0.
         nrml = xyz/norm(xyz)#normalize vector 
         self.nrmls[i] = (info[0],nrml[0],nrml[1],nrml[2])
         nrml = npround(nrml,self.__dp)
         key = (info[0],nrml[0],nrml[1],nrml[2])
         self.rev_nm[key] = i
      return

   def _update_all_nrml(self,vector):
      '''
          update all nrmls with given vector
      '''
      print "vector will be auto normalized"
      assert(len(vector)==3)
      self.rev_nm={}
      for i in self.nrmls:
         info = self.nrmls[i]
         xyz = array(vector,dtype='float64')
         nrml = xyz/norm(xyz)#normalize vector 
         self.nrmls[i] = (info[0],nrml[0],nrml[1],nrml[2])
         nrml = npround(nrml,self.__dp)
         key = (info[0],nrml[0],nrml[1],nrml[2])
         self.rev_nm[key] = i
  


      


   def _rebuild_rev_info(self):
      '''
        assumption: no coincident nodes
        reintialized rev_nd,rev_nm and generate corresponding info based current nodes,nrmls
      '''
      self.rev_nd = {}
      self.rev_nm = {}
      for i in self.nodes:
          xyz = self.nodes[i]
          key = (round(xyz[0],self.__dp),round(xyz[1],self.__dp),round(xyz[2],self.__dp))
          self.rev_nd[key] = i
      for i in self.nrmls:
          info = self.nrmls[i]
          key = (info[0],round(info[1],self.__dp),round(info[2],self.__dp),round(info[3],self.__dp))
          self.rev_nm[key] = i

   def _recreate_avail_info(self):
        '''
           recalcuate __cur_avail info for node/norml/element using current elems/nodes/nrmls
        '''
        self.__cur_avail_el_id = len(self.elems)+1
        self.__cur_avail_nd_id = len(self.nodes)+1
        self.__cur_avail_nm_id = len(self.nrmls)+1


   def _rm_unused_node(self):
        '''
        remove unused node from the model
        '''
        c1 = set()
        l1 =list()
        for i in self.elems:
            tmp = self.elems[i][POS.NODLIST]
            c1.update(tmp)
        for i in self.nodes:
            if not(i in c1):
                l1.append(i)
        for i in l1:
            self.nodes.pop(i)

   def _rm_unused_nrml(self):
        '''
        remove unused nrml from the model
        '''
        c1 = set()
        l1 =list()
        for i in self.elems:
            tmp = self.elems[i][POS.NRMLIST]
            c1.update(tmp)

        for i in self.nrmls:
            if not(i in c1):
                l1.append(i)
        for i in l1:
            self.nrmls.pop(i)

   def _renum_model(self):
        self._rm_unused_node()
        self._rm_unused_nrml()
        # print 'nodes=',len(self.nodes)
        # print 'nrmls=',len(self.nrmls)
        new_nodes={}
        new_nrmls={}
        re_node={}
        re_nrml={}
        ids=1
        ### Assume no coincident nodes/nrmls
        for i in self.nodes:
            new_nodes[ids] = self.nodes[i]
            re_node[i]=ids  
            ids+=1
        
        ids=1
        # print re_node
        for i in self.nrmls:
            info=list(self.nrmls[i])
            info[0]=re_node[info[0]]
            new_nrmls[ids] = tuple(info)
            re_nrml[i]=ids  
            ids+=1

        self.nodes=new_nodes
        self.nrmls=new_nrmls

        for i in self.elems:
            info = self.elems[i]
            nlist=info[POS.NODLIST]
            nrmlist=info[POS.NRMLIST] 
            tmp1=list()
            tmp2=list()
            for j in range(info[POS.TYPE]):
                tmp1.append(re_node[nlist[j]])
                tmp2.append(re_nrml[nrmlist[j]])
            self.elems[i]=[info[0],info[1],tuple(tmp1),tuple(tmp2)]

        self._rebuild_rev_info()
        self._recreate_avail_info()

                




   def check_coincident_nodes(p=False):

        pass

  
   def get_edge_info(self):
      '''
        generate edge information using current elem
      '''
      self.edges = {}
      # self._edge_info = True
      for i in self.elems:
         elem_info = self.elems.get(i,None)
         logging.debug("process elem",i)
         nodelist = elem_info[POS.NODLIST]
         nodelist = list(nodelist)
         nodelist.append(nodelist[0]) 
         logging.debug(nodelist)
         for i in range(elem_info[POS.TYPE]/2):
            e = tuple(sorted([nodelist[2*i],nodelist[2*i+2]]))
            self.edges[e]=nodelist[2*i+1]
       

   def _find_mid_point(self,e):
       ''' 
       if edge e exist in current model, return exisit middle node
       if      e not exist, create new node, update nodes,rev_nd,edge
       ------------------
       return value is ----1) exisiting node id
                           2) new generated node id
       '''
       
       pos = self.edges.get(e,None)
       dp = self.__dp
       if pos == None:
          # NOTE:establish new middle node, but check coincident node
          n1 = e[0]
          n2 = e[1]
          xyz =[0.,0.,0.]
          for i in range(3):
              xyz[i] = (self.nodes[n1][i]+self.nodes[n2][i])/2.
          key = (round(xyz[0],dp),round(xyz[1],dp),round(xyz[2],dp))
          pos = self.rev_nd.get(key,self.__cur_avail_nd_id)
          if pos == self.__cur_avail_nd_id:
              # TODO update rev_nd 
              self.nodes[pos] = tuple(xyz)
              self.rev_nd[key]=pos
              self.edges[e]=pos
              self.__cur_avail_nd_id += 1
       return pos
    
   def horiz_r(self,nid):
       '''
       calculate radius
       --------------
       nid: node id
       '''
       from numpy import sqrt
       n=self.nodes[nid]
       r=sqrt(n[0]**2+n[1]**2)
       return r

   def renew_circular_midpoint(self,ctr=[0.,0.]):
       '''
        regenerate midpoints for circular edge if initally set straight
       '''
       from math import acos
       from numpy import sign
       tol = 10e-7
       for e in self.edges:
           n1 = e[0]
           n2 = e[1]
           if abs(self.nodes[n1][2]-self.nodes[n2][2]) < 1e-7:
                r1 = self.horiz_r(n1)
                r2 = self.horiz_r(n2)
                if (abs(r1-r2)<tol):
                    n3=self._find_mid_point(e)               
                    r3= self.horiz_r(n3)
                    info=list(self.nodes[n3])
                    theta = acos(info[0]/r3)
                    info[0] = sign(info[0])*abs(r1*cos(theta))
                    info[1] = sign(info[1])*abs(r1*sin(theta))
                    self.nodes[n3]=tuple(info)



   ###################################
   # 
   #    Mesh Operation
   #
   ###################################


   def shift_mesh(self,vector):
       '''
       Mesh translation
       '''
       assert(len(vector)==3)
       for i in self.nodes:
           info=list(self.nodes[i])
           for j in range(3):
               info[j]+=vector[j]
           self.nodes[i]=tuple(info)
       self._rebuild_rev_info()
       #FIXME rebuild revnodes

   def scale_mesh(self,factor):
       '''
       Mesh scale
       '''
       assert(factor>0)
       for i in self.nodes:
           info=list(self.nodes[i])
           for j in range(3):
               info[j]=info[j]*factor
           self.nodes[i]=tuple(info) 
       self._rebuild_rev_info

   # def mirror_mesh(self,kind=2,base=[0,0,0]):
   # TODO there is direction problem with this function
       # """mirror_mesh(kind=2,base=[0,0,0])
         # kind=0,along yz plane
         # kind=1,along xz plane
         # kind=2,along xz plane
       # """
       # n = deepcopy(self)
       # for i in n.nodes:
           # info = list(n.nodes[i])
           # info[kind] = 2*base[kind]-info[kind] 
           # n.nodes[i]=tuple(info)
       # for i in n.nrmls:
           # info = list(n.nrmls[i])
           # info[kind+1] =-info[kind+1] 
           # n.nrmls[i]=tuple(info)
       # n._rebuild_rev_info()
       # return n
   
   # def combine_mesh(m1,m2):
      # '''
      # n = m1.combine_mesh(m1,m2)
      # combine m1,m2 and return the new mesh
      # '''
      # m3 = deepcopy(m1)
      # m3.devour_mesh(m2)
      # return m3

   def devour_mesh(self, new_mesh):
      '''
       n.devour_mesh(m1)
       add m1 to exisiting mesh n
      '''
      dp = self.__dp
      renum_node ={}
      for i in new_mesh.nodes:
         xyz = new_mesh.nodes[i]
         key = (round(xyz[0],dp),round(xyz[1],dp),round(xyz[2],dp))
         pos = self.rev_nd.get(key,self.__cur_avail_nd_id)
         renum_node[i]=pos
         if pos == self.__cur_avail_nd_id:
            # create new id, record in rev_nd
            self.nodes[pos] = xyz
            self.rev_nd[key] = pos
            self.__cur_avail_nd_id+=1
      logging.debug(renum_node)
      renum_nrml={}
      for i in new_mesh.nrmls:
         nlist = list(new_mesh.nrmls[i])
         # logging.debug(nlist)
         nlist[0] = renum_node[nlist[0]]
         # logging.debug(nlist)
         key = tuple(nlist)
         pos = self.rev_nm.get(key,self.__cur_avail_nm_id)
         renum_nrml[i] = pos
         if pos == self.__cur_avail_nm_id:
            self.nrmls[pos] = key
            self.rev_nm[key] = pos
            self.__cur_avail_nm_id+=1
      # logging.debug(renum_nrml) 
      # NOTE:no element coincidence check
      for i in new_mesh.elems:
         einfo = new_mesh.elems[i]
         nodelist = list(einfo[POS.NODLIST])
         nrmlist = list(einfo[POS.NRMLIST])
         # logging.debug(nrmlist)
         for j in range(einfo[POS.TYPE]):
            nodelist[j] = renum_node[nodelist[j]]
            nrmlist[j] = renum_nrml[nrmlist[j]]
         # logging.debug(nrmlist)
         self.elems[self.__cur_avail_el_id] = [einfo[0],einfo[1],tuple(nodelist),tuple(nrmlist)]
         self.__cur_avail_el_id+=1
            
         

   def __test(self):
        print "i am here"
        pass



    ###################################################################
    #                                                                 #
    #                INPUT AND OUTPUT PART OF CLASS                   #
    #                                                                 #
    ###################################################################
   def output_mesh(self,path,kind):
       """
       output_mesh(path,kind)
       kind = 0, body mesh
       kind = 1, no damp info
       kind = 2, damp info =0,mesh do not need to have damp info setup
       kind = 3, need damp info setup
       currently: only 8-node elem
       """
       if kind == FILETYPE.BODY:
           self.__output_as_bd(path)
       if kind == FILETYPE.SURFACE_WO_DAMP:
           self.__output_as_fs(path,0)
       if kind == FILETYPE.SURFACE_W_DAMP:
           self.__output_as_fs(path,2)
       if kind == FILETYPE.SURFACE_AUTO_DAMP:
           self.__output_as_fs(path,1)
           
   def __output_as_fs(self,path,kind=1):
       """
       _output_as_fs(path,kind=1)
       kind = 0, no damp info
       kind = 1, damp info =0,mesh do not need to have damp info setup
       kind = 2, need damp info setup
       currently: only 8-node elem
       """
       print "only 8-node elem has been implemented"
       
       f_sf = open(path,"w")
       f_sf.write('{0:<5d}'.format(self.__cur_avail_el_id-1) \
               + '{0:<5d}\n'.format(self.__cur_avail_nd_id-1))
       acc_sf_elem = 1
       for i_elem in self.elems:
            str1 = ''
            str2 = ''
            str3 = ''
            nodelist = self.elems[i_elem][POS.NODLIST]
            for j in range(8):#TODO add triangle support
                    
                    str1 += '{0:<9.6f}     '.format(self.nodes[nodelist[j]][0])
                    str2 += '{0:<9.6f}     '.format(self.nodes[nodelist[j]][1])
                    if kind==2:
                        #NOTEME damp info is not defined in this class
                        assert(hasattr(self,'damp_info'))
                        dval = self.damp_info.get(nodelist[j],0.)
                        str3 += '{0:<9.6f}     '.format(dval)
                    if kind==1:
                        str3 += '{0:<9.6f}     '.format(0.)

            f_sf.write(('{0:<6d}   8\n').format(acc_sf_elem))
            f_sf.write(str1+'\n')
            f_sf.write(str2+'\n')
            if kind:
                f_sf.write(str3+'\n')
            acc_sf_elem += 1
       f_sf.close()

   def __output_as_bd(self,path):
       """
       _output_as_bd(path)
       output using body mesh format
       """
       f = open(path,"w")
       f.write('0\n') # write isys
       # READ(2,*)   NELEMB, NNB, NNBD, IPOL
       f.write('    '.join([' ',str((len(self.elems))),str(len(self.nodes)),str(len(self.nrmls)),'1']))
       f.write("\n")
       f.write('1  0  0.00  0.00  0.00')# FIXME second maybe 0
       f.write("\n")

       for i_node in self.nodes:
               f.write(('{0:<7d}'.format(i_node)))
               f.write('1     ')
               f.write('    '.join('{0:<9.4f}'.format(i) for i in self.nodes[i_node]))
               f.write("\n")

       for i_norm in self.nrmls:
               f.write(('{0:<7d}'.format(i_norm)))
               f.write('1     ')
               f.write('    '.join('{0:<9.4f}'.format(i) for i in self.nrmls[i_norm][1:4]))
               f.write("\n")

       for i_elem in self.elems:
               f.write(('{0:<7d}   8\n').format(i_elem))
               #f.write(('{0:<5d}'.format(i_elem)))
               f.write('    '.join(str(i) for i in self.elems[i_elem][POS.NODLIST]))
               f.write("\n")

       for i_elem in self.elems:
               f.write(('{0:<7d}   8\n').format(i_elem))
               #f.write(('{0:<5d}'.format(i_elem)))
               f.write('    '.join(str(i) for i in self.elems[i_elem][POS.NRMLIST]))
               f.write("\n")

       f.close()


    ###################################################################
    #                                                                 #
    #                INPUT AND OUTPUT PART OF CLASS                   #
    #                                                                 #
    ###################################################################
   def read_mesh(self,path,kind,vector=[0,0,0]):
       '''
           def read_mesh(self,path,kind,vector=[0,0,0])
           
           kind=0 for body
           kind=1 for surface w/o damp info
           kind=2 for surface w   damp info
           kind=9 for external 4node data
           vector apply to external only

       '''
       if kind == FILETYPE.BODY:
           self.__read_body_fmt(path)
       if kind == FILETYPE.SURFACE_WO_DAMP:
           self.__read_surface_fmt(path,False)
       if kind == FILETYPE.SURFACE_W_DAMP:
           self.__read_surface_fmt(path,True)
       if kind == FILETYPE.EXTERNAL:
           self.__read_external(path,[0,0,0]) 


   def __read_surface_fmt(self,path,flag_damp=True):
       '''
          Suppose surface only have 8 nodes element
       '''

       with open(path,"r") as f:
          num_elem =[int(i) for i in f.readline().split()][0]
          dp = self.__dp
          for i_elem in range(num_elem):
              tmp = f.readline().split()
              tmp1 = [float(i) for i in f.readline().split()]
              tmp2 = [float(i) for i in f.readline().split()]
              if flag_damp==True:
                  tmp3 = [float(i) for i in f.readline().split()]

              nodelist = []
              for j in range(8):
                  node = (round(tmp1[j],dp),round(tmp2[j],dp),round(0.0,dp))
                  if node in self.rev_nd:
                      pos = self.rev_nd[node]
                      nodelist.append(pos)
                  else:
                      self.rev_nd[node] = self.__cur_avail_nd_id#create new node in rev dict
                      self.nodes[self.__cur_avail_nd_id] = (tmp1[j],tmp2[j],0.0)#create new node in dict_fs_node
                      nodelist.append(self.__cur_avail_nd_id)
                      if flag_damp:
                          self.damp_info[self.__cur_avail_nd_id] = tmp3[j]
                      self.__cur_avail_nd_id+=1
              self.elems[self.__cur_avail_el_id] = ['free surface',8,tuple(nodelist),tuple(nodelist)]
              self.__cur_avail_el_id += 1
          for i_node in self.nodes:
              key = (i_node,round(0.0,self.__dp),round(0.0,self.__dp),round(1.,self.__dp))
              self.nrmls[self.__cur_avail_nm_id] = key 
              self.rev_nm[key] = self.__cur_avail_nm_id
              self.__cur_avail_nm_id+=1
            


   def __read_body_fmt(self,path):
        dp = self.__dp
        f = open(path,'r')
        flag =[int(i) for i in f.readline().split()][0]
        # read number of elem,node,normal, axis
        tmp = [int(i) for i in f.readline().split()]
        num_bd_elem = tmp[0]
        num_bd_node = tmp[1]
        num_bd_nrml = tmp[2]
        num_bd_axis = tmp[3]
        ####----------header finished-----
        dict_axis = {}
        for i in range(num_bd_axis):
            tmp = [float(i) for i in f.readline().split()]
            dict_axis[int(tmp[0])] = (int(tmp[1]),tmp[2:5]) 
            # first number: index ,sec : 1: polar,0:cartesian ,third-fifth: offset coordinates
        ##--------coordinate info finished----------

        renum_bd_node = {}#key is newly assigned id, value is old id
        for i_node in range(num_bd_node): 
            tmp = [float(i) for i in f.readline().split()]
            offset = dict_axis[int(tmp[1])][1]
            axis_type = dict_axis[int(tmp[1])][0]
            if(axis_type==0):
                x = tmp[2]+offset[0]
                y = tmp[3]+offset[1]
            if(axis_type==1):
                x = tmp[2]*cos(deg2rad(tmp[3]))+offset[0]
                y = tmp[2]*sin(deg2rad(tmp[3]))+offset[1]
            key = (round(x,dp),round(y,dp),round(tmp[4],dp))
            
            if (key in self.rev_nd):
                pos = self.rev_nd[key]
                renum_bd_node[int(tmp[0])] = pos         

            else:
                self.rev_nd[key] = self.__cur_avail_nd_id
                renum_bd_node[int(tmp[0])] = self.__cur_avail_nd_id         
                self.nodes[self.__cur_avail_nd_id] = (x,y,tmp[4])
                self.__cur_avail_nd_id+=1


        ####-----------------------------------------
        # processing normal data    
        '''
           read in nrml data, saved in tmp database
           no need renumbering, no need create loop-up dic
        '''
        tmp_nrmls={}
        for i_nrml in range(num_bd_nrml):
            tmp = [float(i) for i in f.readline().split()]
            offset = dict_axis[int(tmp[1])][1]
            axis_type = dict_axis[int(tmp[1])][0]
            if(axis_type==0):
                x = tmp[2]+offset[0]
                y = tmp[3]+offset[1]
            if(axis_type==1):
                x = tmp[2]*cos(deg2rad(tmp[3]))+offset[0]
                y = tmp[2]*sin(deg2rad(tmp[3]))+offset[1]
            ######--------------------
            tmp_nrmls[int(tmp[0])] = (x,y,tmp[4])   

        ###############NOTE:read in elem info---
        renum_elem={}
        for i_elem in range(num_bd_elem):
            tmp = [int(i) for i in f.readline().split()]
            tmp1 = [int(i) for i in f.readline().split()]
            nodelist = []
            renum_elem[tmp[0]]=self.__cur_avail_el_id
            ## update node number in nodelist
            for k in range(tmp[1]):
                nodelist.append(renum_bd_node[tmp1[k]])
            self.elems[self.__cur_avail_el_id] = ['body',tmp[1],tuple(nodelist)]   
            nodelist.append(nodelist[0])#for use(n,n+1) edge pair
            self.__cur_avail_el_id+=1

        #===============================================================
        tmp_elem_nrml={}# tmporialy save nrmlist info
        for i_elem in range(num_bd_elem):
            tmp = [int(i) for i in f.readline().split()]
            tmp1 = [int(i) for i in f.readline().split()]
            tmp_elem_nrml[renum_elem[i_elem+1]]=tmp1[0:9]  # pair nrmlist with renumber elem id 
            
        '''
           for each element, associate node with corresponding nrml
           renumber nrml, find associated note create self.nrmls dict,update
           self.elems with nrmlist
        '''
        renum_nrml = {}
        processed_nrml=[]
        for ie in tmp_elem_nrml:
            nrmlist = tmp_elem_nrml[ie]
            for j in range(len(nrmlist)):
                if nrmlist[j] in processed_nrml:
                   nrmlist[j] = renum_nrml[nrmlist[j]]# update id in old nrmlist
                else:
                   anode = self.elems[ie][2][j] #associated node = jth node in node list 
                   cnrml = nrmlist[j]
                   pos = tmp_nrmls[cnrml]
                   key = (anode,round(pos[0],self.__dp),round(pos[1],self.__dp),round(pos[2],self.__dp))
                   
                   loc = self.rev_nm.get(key,self.__cur_avail_nm_id)
                   renum_nrml[cnrml] = loc# record id change
                   #FIXME NO there one to multiple case
                   nrmlist[j] = loc# update id in old nrmlist
                   processed_nrml.append(cnrml)
                   
                   if loc == self.__cur_avail_nm_id: 
                      self.nrmls[self.__cur_avail_nm_id]=key#record new nrml
                      self.rev_nm[key] = self.__cur_avail_nm_id#record new key
                      self.__cur_avail_nm_id+=1
                   else:
                      pass
            self.elems[ie].append(tuple(nrmlist))#append updated nrmlist

   def __read_external(self,path,vector,z_offset=0.,tag='external'):
       print "pls make sure edge info is ready,if you want to use old middle points"
       assert(isinstance(vector,list))
       assert(len(vector)==3)
      # vector is nrml vector,kind is elem kind('free','body',or user defined)
      
       dp = self.__dp
   #     num_add_node = 12
       # num_add_elem = 6

       node_set=set()
      
       #TODO complete this, make sure edge_info is generated before read add
       f = open(path,"r")
       renum_add_node = {}
       # for i_node in range(num_add_node): 
   
       while True:
           tmp = [float(i) for i in f.readline().replace('elem','0').split()]
           if len(tmp) != 4:
               break
           x = tmp[1]
           y = tmp[2]
           z = tmp[3]+z_offset    
           key = (round(x,dp),round(y,dp),round(z,dp))
           pos = self.rev_nd.get(key,self.__cur_avail_nd_id)
           renum_add_node[int(tmp[0])] = pos         
           node_set.add(pos)
           if pos == self.__cur_avail_nd_id:
               self.rev_nd[key] = pos
               self.nodes[pos] = (x,y,z)
               self.__cur_avail_nd_id+=1
       logging.debug(renum_add_node)
       logging.debug(len(renum_add_node))

       flag = True
       tmp = [int(i) for i in tmp]#first line for elem already read
       renum_elem={}

       # TODO throw duplicate elems
       # for i_elem in range(num_add_elem):
       while flag == True:
           print tmp
           nodelist = []
           ## update node number in nodelist
           for k in range(4):#FIXME add traingle support
               nodelist.append(renum_add_node[tmp[k+2]])
               nodelist.append(0)
           nodelist.append(nodelist[0])
           for e in range(4):#FIXME add triangle support
               edge = tuple(sorted([nodelist[2*e],nodelist[2*e+2]]))
               midp = self._find_mid_point(edge)
               nodelist[2*e+1] = midp
               node_set.add(midp)

               # renum_add_node[midp]=midp#NOTE {not work with coincident node for input}this fixed nrml generation problem
           nodelist.pop(8)
           # logging.debug(nodelist)
           # nrmlist = np.array(nodelist)+self.__cur_avail_nm_id-1# nrmlist id shoud be offset according to current num of nrmls
           renum_elem[tmp[0]]=self.__cur_avail_el_id
           self.elems[self.__cur_avail_el_id] = [tag,8,tuple(nodelist)]#,tuple(nrmlist)]   
           self.__cur_avail_el_id+=1
           tmp = [int(i) for i in (f.readline().replace('elem','0')).split()]
           if len(tmp)==0:
               flag = False



       
       node_2_nrml={}
       for i_node in node_set :
           key=(i_node,round(vector[0],dp),round(vector[1],dp),round(vector[2],dp))
           # NOTE didn't check coincident nrml
           self.nrmls[self.__cur_avail_nm_id] = key
           self.rev_nm[key] = self.__cur_avail_nm_id
           node_2_nrml[i_node] = self.__cur_avail_nm_id

           self.__cur_avail_nm_id+=1
       # logging.debug(renum_elem)
       logging.debug(node_2_nrml)
       for i_elem in renum_elem:
           elem = self.elems[renum_elem[i_elem]]
           nodelist = elem[2]
           nrmlist=[]
           logging.debug(nodelist)
           for j in range(elem[1]):
               logging.debug(j)
               nrmlist.append(node_2_nrml[nodelist[j]])
           self.elems[renum_elem[i_elem]].append(nrmlist)


       
      #TODO update nrmls and rev_nrmls
      

