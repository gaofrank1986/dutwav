import numpy as np
class FreeTerm(object):
    def __init__(self):
        self.result = {}
        # self.result2 ={}
        self.__dp = 4
        self.list1=[]#frc31 > threshold
        self.list2=[]#frc32 > threshold
        self.info={}
        # self.info2={}

    def read_result(self,path,dp=4):
        list1=[]
        list2=[]
        with open(path,"rb") as f:
            lines = f.readlines()
            acc = 0
            for line in lines:
                acc+=1
                tmp = [float(i) for i in line.split()]
                # key = tuple(np.round(tmp[0:2],dp))
                key = (round(tmp[0],dp),round(tmp[1],dp))
                self.result[key] = tmp[2:]
                self.info[acc] = key 
                # self.info2[key] = acc
                # self.result2[acc] = tmp[2:6]
                if abs(tmp[3])>1e-2:
                    list1.append(acc)
                if abs(tmp[4])>1e-2:
                    list2.append(acc)
        self.list1=list1
        self.list2=list2
    
    def calc_fterm(self,r2):
        for (pos,fterm) in self.result.items():
            if pos in r2.result.keys():
                f2 = r2.result[pos]
                fterm[0] = 0. - f2[0]
                fterm[1] = 0. - f2[1]
                fterm[2] = 0. - f2[2]
                fterm[3] = f2[3]# NOTE 1- is taken away

            else:
                fterm[0] = 0.
                fterm[1] = 0.
                fterm[2] = 0.
                fterm[3] = 0.5
            self.result[pos][0:4]=fterm[0:4]

    def print_fterm(self):
        for (i,key) in self.info.items():
            print self.result[key]

    def output_fterm(self,path):
        with open(path,"wb") as f:
            for (i,key) in self.info.items():
                f.write(('{0:<5d}'.format(i)))
                f.write('    '.join('{0:<12.8f}'.format(j) for j in self.result[key]))
                f.write("\n")



class AnalyticalPotential(object):
    def __init__(self):
        self.vec=None
        self.eti=None

    def read_vector(self,path):
        with open(path,"rb") as f:
            tmp = f.readlines()
            vec = np.zeros((len(tmp),1),dtype = 'complex128')
            for i in range(len(tmp)):
                s = tmp[i]
                # s = s.replace('(',' ')
                # s = s.replace(')',' ')
                # s2=[float(j) for j in s.split(',')]
                s2=[float(j) for j in s.split()]
                vec[i,0] = s2[1]+1j*s2[2]
            self.vec=vec
    def get_eti(self,omega):
        # eti = iw/g * potential
        self.eti = self.vec*1j*omega/9.807

    def output(self,path):
        re= np.zeros((self.vec.shape[0],2),dtype='float64')
        re[:,0]=self.vec.real[:,0]
        re[:,1]=self.eti.real[:,0]
        np.savetxt(path,re)



class BoundaryValue(object):
    def __init__(self):
        self.bc=None

    def read_bc(self,path):
        with open(path,"rb") as f:
            tmp = f.readlines()
            vec = np.zeros((len(tmp),1),dtype = 'float64')
            for i in range(len(tmp)):
                vec[i,0] = float(tmp[i])
            self.bc=vec

class MatrixA(object):
    def __init__(self):
        self.C=None
        self.A=None
        self.n=0
        self.m=0
        self.err={}
        self.derr={}
    def read_matrix(self,path):
        self.A = np.loadtxt(path,dtype='float64')
        n = int(np.sqrt(self.A.size))
        self.A.resize(n,n)
        self.n=n

    def read_matrix_c(self,path):
        self.C = np.loadtxt(path,dtype='float64')
        m = int(self.C.size/self.n)
        self.C.resize(self.n,m)
        self.m=m

    def analysis_matrix(self,threshold):
        import copy
        self.err={}
        self.derr={}
        b = copy.copy(self.A)
        while (b.max()) >threshold:
            n1 = b.argmax()/self.n
            n2 = b[n1,:].argmax()
            key = (n1+1,n2+1)
            if n1!=n2:
                self.err[key]=b[n1,n2]
            else:
                self.derr[key]=b[n1,n2]
            b[n1,n2]=0.


class ElemValue():
    def __init__(self):
        self.res={}

    def read_result(self,path):
        with open(path,"rb") as f:
            acc=0
            for line in f:
                line = line.split()
                acc+=1
                value=[float(j) for j in line[3:]]
                eid = int(line[1])
                nid = int(line[0])
                self.res[nid]=self.res.get(nid,{})
                self.res[nid][eid]=self.res[nid].get(eid,{})
                flag = int(line[2])
                if acc==1:
                    self.res[nid][eid]['amat']=value
                if acc==2:
                    self.res[nid][eid]['bmat']=value
                if acc==3:
                    self.res[nid][eid]['aval']=value
                if acc==4:
                    self.res[nid][eid]['bval']=value
                    acc=0
    def pp_print(self,nid,eid):
        value = None
        try:
            value = self.res[nid][eid]
        except:
            print "Error Accessing, Check input number"
        if value:
            for i in value:
                str1='{0:5d} {0:5d}    '.format(nid,eid)
                str1+=i
                str1+="   "
                for j in range(8):
                    str1 += "{0:14.8f}".format(value[i][j])
                print str1





