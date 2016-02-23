import logging
import numpy as np
from enum import Enum


class POS(Enum):
    TAG=0
    TYPE=1
    NODLIST=2
    NRMLIST=3

class Node(object):
    # __slots__={'pos','id'}
    def __init__(self,id,pos):
        self.pos= pos
        self.id = id
        self.connected_elem=[]
        self.connected_edge=[]
        self.associated_nrml=[]
        self.coor = None
    
    def print_node(self):
        pass

class Normal(object):
    # __slots__={'pos','id'}
    def __init__(self,id,pos):
        self.pos= pos
        self.id = id
        self.associated_node=[]
    
    def print_nrml(self):
        pass

class Coordinate(object):
    def __init__(self,id,ctype,offset):
        self.id = id
        self.type = ctype
        self.offset = offset

    def print_coor():
        pass

# ====================================
class Element(object):
    # __slots__={'nl','nrl','etype','kind','id'}
   def __init__(self,ename):
      self.tag = tag
      self.etype = -1#6 for triangle,8 for rectangular,4 for 4 node rect,3 for 3 node triangle
      self.nl=[]
      self.nrl=[]

   def print_elem(self,flag=0): # flag =0  for nodelist,flag=1 for nrml list
        pass

class _Mesh_core(object):
   def __init__(self,dp=4):
        self.nsys = 0
        self.nelem = 0
        self.nnode = 0
        self.nnrml = 0
        

        self.nodes = {}
        self.elems = {}
        self.rev_nd = {}
        self.rev_nm = {}
        self.rev_el={}
        self.nrmls ={}
        self.coors = {}
        self.edge_mdp={}
        # self.damp={}
        
        self._cur_avail_nd_id=1
        self._cur_avail_el_id=1
        self._cur_avail_nm_id=1
        self._dp=4

   def is_exist_node(self,pos):
        assert(len(pos)==3)
        key = (round(pos[0],self._dp),round(pos[1],self._dp),round(pos[2],self._dp))
        return(key in self.rev_nd)

   def is_exist_nrml(self,node_id,pos):
       assert(len(pos)==3)
       key = (node_id,round(pos[0],self._dp),round(pos[1],self._dp),round(pos[2],self._dp))
       return(key in self.rev_nd)



   def _is_node_fs(self,n,z=0):
       xyz=self.nodes[n]
       if abs(xyz[2]-z)<1e-4:
           return True
       return False
   
   def test_mesh(self):
      """
      first n node must be surface node
      first elems must be surface elems
      """
      pass

        
#    def _count_elem(self,tag=[]):
       # result={}
       # for i in self.elems:
           # key = self.elems[i][0]
           # result[key] = result.get(key,0)+1
       # return result 
           
   def get_edge_info(self):
      self.edge_mdp = {}
      self._edge_info = True
      for i in self.elems:
         elem_info = self.elems.get(i,None)
         logging.debug("process elem",i)
         nodelist = elem_info[POS.NODLIST]
         nodelist = list(nodelist)
         nodelist.append(nodelist[0]) 
         logging.debug(nodelist)
         for i in range(elem_info[POS.TYPE]/2):
            e = tuple(sorted([nodelist[2*i],nodelist[2*i+2]]))
            self.edge_mdp[e]=nodelist[2*i+1]
       


#    def get_waterline_node(self):
       # print "elements in mesh must have well-defined tag\n\
       # only 'body' elements will be looked into"
       # nodes=set()
       # for i in self.elems:
           # info = self.elems[i]
           # if info[0]=='body':
               # nodelist = info[POS_NODLIST]
               # for j in range(info[POS_ELEMTYPE]):
                   # if self._is_node_fs(nodelist[j]):
                       # nodes.add(nodelist[j])
       # return(nodes)                
                    
   def _rebuild_rev_info(self):
      self.rev_nd = {}
      self.rev_nm = {}
      for i in self.nodes:
          xyz = self.nodes[i]
          key = (round(xyz[0],self._dp),round(xyz[1],self._dp),round(xyz[2],self._dp))
          self.rev_nd[key] = i
      for i in self.nrmls:
          info = self.nrmls[i]
          key = (info[0],round(info[1],self._dp),round(info[2],self._dp),round(info[3],self._dp))
          self.rev_nm[key] = i
       


   def combine_mesh(m1,m2):
      '''
      n = m1.combine_mesh(m1,m2)
      combine m1,m2 and return the new mesh
      '''
      import copy
      m3 = copy.deepcopy(m1)
      m3.devour_mesh(m2)
      return m3

   def devour_mesh(self, new_mesh):
      '''
       n.devour_mesh(m1)
       add m1 to exisiting mesh n
      '''
      dp = self._dp
      renum_node ={}
      for i in new_mesh.nodes:
         xyz = new_mesh.nodes[i]
         key = (round(xyz[0],dp),round(xyz[1],dp),round(xyz[2],dp))
         pos = self.rev_nd.get(key,self._cur_avail_nd_id)
         renum_node[i]=pos
         if pos == self._cur_avail_nd_id:
            # create new id, record in rev_nd
            self.nodes[pos] = xyz
            self.rev_nd[key] = pos
            self._cur_avail_nd_id+=1
      logging.debug(renum_node)
      renum_nrml={}
      for i in new_mesh.nrmls:
         nlist = list(new_mesh.nrmls[i])
         # logging.debug(nlist)
         nlist[0] = renum_node[nlist[0]]
         # logging.debug(nlist)
         key = tuple(nlist)
         pos = self.rev_nm.get(key,self._cur_avail_nm_id)
         renum_nrml[i] = pos
         if pos == self._cur_avail_nm_id:
            self.nrmls[pos] = key
            self.rev_nm[key] = pos
            self._cur_avail_nm_id+=1
      # logging.debug(renum_nrml) 
      # NOTE:no element coincidence check
      for i in new_mesh.elems:
         einfo = new_mesh.elems[i]
         nodelist = list(einfo[POS.NODELIST])
         nrmlist = list(einfo[POS.NODLIST])
         # logging.debug(nrmlist)
         for j in range(einfo[POS.TYPE]):
            nodelist[j] = renum_node[nodelist[j]]
            nrmlist[j] = renum_nrml[nrmlist[j]]
         # logging.debug(nrmlist)
         self.elems[self._cur_avail_el_id] = (einfo[0],einfo[1],tuple(nodelist),tuple(nrmlist))
         self._cur_avail_el_id+=1
            
         
         



   def output_as_tecplt(self,path):
       pass

               

       


   def _find_mid_point(self,e):
       ''' 
       if edge e exist in current model, return exisit middle node
       if      e not exist, create new node, update nodes,rev_nd,edge
       '''
       
       pos = self.edge_mdp.get(e,None)
       dp = self._dp
       if pos == None:
          # NOTE:establish new middle node, but check coincident node
          n1 = e[0]
          n2 = e[1]
          xyz =[0.,0.,0.]
          for i in range(3):
              xyz[i] = (self.nodes[n1][i]+self.nodes[n2][i])/2.
          key = (round(xyz[0],dp),round(xyz[1],dp),round(xyz[2],dp))
          pos = self.rev_nd.get(key,self._cur_avail_nd_id)
          if pos == self._cur_avail_nd_id:
              # TODO update rev_nd 
              self.nodes[pos] = tuple(xyz)
              self.rev_nd[key]=pos
              self.edge_mdp[e]=pos
              self._cur_avail_nd_id += 1
       return pos



    ###################################################################
    #                                                                 #
    #                INPUT AND OUTPUT PART OF CLASS                   #
    #                                                                 #
    ###################################################################
         
   def _output_as_fs(self,path,kind=1):
       """
       _output_as_fs(path,kind=1)
       kind = 0, no damp info
       kind = 1, damp info =0,mesh do not need to have damp info setup
       kind = 2, need damp info setup
       currently: only 8-node elem
       """
       f_sf = open(path,"w")
       f_sf.write('{0:<5d}'.format(self._cur_avail_el_id-1) + '{0:<5d}\n'.format(self._cur_avail_nd_id-1))
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
                        dval = self.damp_info.get(nodelist[j],0.)
                        str3 += '{0:<9.6f}     '.format(dval)
                    if kind==1:
                        str3 += '{0:<9.6f}     '.format(0.)

            f_sf.write(('{0:<5d}   8\n').format(acc_sf_elem))
            f_sf.write(str1+'\n')
            f_sf.write(str2+'\n')
            if kind:
                f_sf.write(str3+'\n')
            acc_sf_elem += 1
       f_sf.close()

   def _output_as_bd(self,path):
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
               f.write(('{0:<5d}'.format(i_node)))
               f.write('1     ')
               f.write('    '.join('{0:<9.4f}'.format(i) for i in self.nodes[i_node]))
               f.write("\n")

       for i_norm in self.nrmls:
               f.write(('{0:<5d}'.format(i_norm)))
               f.write('1     ')
               f.write('    '.join('{0:<9.4f}'.format(i) for i in self.nrmls[i_norm][1:4]))
               f.write("\n")

       for i_elem in self.elems:
               f.write(('{0:<5d}   8\n').format(i_elem))
               #f.write(('{0:<5d}'.format(i_elem)))
               f.write('    '.join(str(i) for i in self.elems[i_elem][POS_NODLIST]))
               f.write("\n")

       for i_elem in self.elems:
               f.write(('{0:<5d}   8\n').format(i_elem))
               #f.write(('{0:<5d}'.format(i_elem)))
               f.write('    '.join(str(i) for i in self.elems[i_elem][POS_NRMLIST]))
               f.write("\n")

       f.close()


    ###################################################################
    #                                                                 #
    #                INPUT AND OUTPUT PART OF CLASS                   #
    #                                                                 #
    ###################################################################

   def read_sf_file(self,path,flag_damp=True):
       print "pls set flag_damp=True if has damp info"

       with open(path,"r") as f:
          num_elem =[int(i) for i in f.readline().split()][0]
          dp = self._dp
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
                      self.rev_nd[node] = self._cur_avail_nd_id#create new node in rev dict
                      self.nodes[self._cur_avail_nd_id] = (tmp1[j],tmp2[j],0.0)#create new node in dict_fs_node
                      nodelist.append(self._cur_avail_nd_id)
                      if flag_damp:
                          self.damp_info[self._cur_avail_nd_id] = tmp3[j]
                      self._cur_avail_nd_id+=1
              self.elems[self._cur_avail_el_id] = ['free surface',8,tuple(nodelist),tuple(nodelist)]
              self._cur_avail_el_id += 1
          for i_node in self.nodes:
              key = (i_node,round(0.0,self._dp),round(0.0,self._dp),round(1.,self._dp))
              self.nrmls[self._cur_avail_nm_id] = key 
              self.rev_nm[key] = self._cur_avail_nm_id
              self._cur_avail_nm_id+=1
              #TODO assign num of nrmls,elems,nodes,zone_list
            


   def read_bd_file(self,path):
        from numpy import cos,sin
        dp = self._dp
        f = open(path,'r')
        flag =[int(i) for i in f.readline().split()][0]
        #TODO  assign self. nsys
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
        # node_waterline = []
        for i_node in range(num_bd_node): 
            tmp = [float(i) for i in f.readline().split()]
            offset = dict_axis[int(tmp[1])][1]
            axis_type = dict_axis[int(tmp[1])][0]
            if(axis_type==0):
                x = tmp[2]+offset[0]
                y = tmp[3]+offset[1]
            if(axis_type==1):
                x = tmp[2]*cos(np.deg2rad(tmp[3]))+offset[0]
                y = tmp[2]*sin(np.deg2rad(tmp[3]))+offset[1]
            key = (round(x,dp),round(y,dp),round(tmp[4],dp))
            
            if (key in self.rev_nd):
                pos = self.rev_nd[key]
                renum_bd_node[int(tmp[0])] = pos         

                # if (abs(tmp[4])<1e-4):
                    # node_waterline.append(pos)

            else:
                self.rev_nd[key] = self._cur_avail_nd_id
                renum_bd_node[int(tmp[0])] = self._cur_avail_nd_id         
                self.nodes[self._cur_avail_nd_id] = (x,y,tmp[4])
                #####NOTE: get waterline node id list
                # if (abs(tmp[4])<1e-4):
                    # node_waterline.append(self._cur_avail_nd_id)
                #### finish water line node
                self._cur_avail_nd_id+=1

        # self._waterline = list(set(node_waterline))


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
                x = tmp[2]*cos(np.deg2rad(tmp[3]))+offset[0]
                y = tmp[2]*sin(np.deg2rad(tmp[3]))+offset[1]
            ######--------------------
            tmp_nrmls[int(tmp[0])] = (x,y,tmp[4])   

        ###############NOTE:read in elem info---
        renum_elem={}
        for i_elem in range(num_bd_elem):
            tmp = [int(i) for i in f.readline().split()]
            tmp1 = [int(i) for i in f.readline().split()]
            nodelist = []
            renum_elem[tmp[0]]=self._cur_avail_el_id
            ## update node number in nodelist
            for k in range(tmp[1]):
                nodelist.append(renum_bd_node[tmp1[k]])
            self.elems[self._cur_avail_el_id] = ['body',tmp[1],tuple(nodelist)]   
            nodelist.append(nodelist[0])#for use(n,n+1) edge pair
            self._cur_avail_el_id+=1

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
                   key = (anode,round(pos[0],self._dp),round(pos[1],self._dp),round(pos[2],self._dp))
                   
                   loc = self.rev_nm.get(key,self._cur_avail_nm_id)
                   renum_nrml[cnrml] = loc# record id change
                   #FIXME NO there one to multiple case
                   nrmlist[j] = loc# update id in old nrmlist
                   processed_nrml.append(cnrml)
                   
                   if loc == self._cur_avail_nm_id: 
                      self.nrmls[self._cur_avail_nm_id]=key#record new nrml
                      self.rev_nm[key] = self._cur_avail_nm_id#record new key
                      self._cur_avail_nm_id+=1
                   else:
                      pass
            self.elems[ie].append(tuple(nrmlist))#append updated nrmlist

   def read_added(self,path,vector,z_offset=0.,comment='add mesh'):
      # vector is nrml vector,kind is elem kind('free','body',or user defined)
      
       dp = self._dp
       num_add_node = 12
       num_add_elem = 6

       node_set=set()
      
       assert(len(vector)==3)
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
           pos = self.rev_nd.get(key,self._cur_avail_nd_id)
           renum_add_node[int(tmp[0])] = pos         
           node_set.add(pos)
           if pos == self._cur_avail_nd_id:
               self.rev_nd[key] = pos
               self.nodes[pos] = (x,y,z)
               self._cur_avail_nd_id+=1
       logging.debug(renum_add_node)
       logging.debug(len(renum_add_node))

       flag = True
       tmp = [int(i) for i in tmp]
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
           # nrmlist = np.array(nodelist)+self._cur_avail_nm_id-1# nrmlist id shoud be offset according to current num of nrmls
           renum_elem[tmp[0]]=self._cur_avail_el_id
           self.elems[self._cur_avail_el_id] = [comment,8,tuple(nodelist)]#,tuple(nrmlist)]   
           self._cur_avail_el_id+=1
           tmp = [int(i) for i in (f.readline().replace('elem','0')).split()]
           if len(tmp)==0:
               flag = False
       
       node_2_nrml={}
       for i_node in node_set :
           key=(i_node,round(vector[0],dp),round(vector[1],dp),round(vector[2],dp))
           # NOTE didn't check coincident nrml
           self.nrmls[self._cur_avail_nm_id] = key
           self.rev_nm[key] = self._cur_avail_nm_id
           node_2_nrml[i_node] = self._cur_avail_nm_id

           self._cur_avail_nm_id+=1
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
      

