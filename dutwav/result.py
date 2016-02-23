class freeterm_rslt(object):
    def __init__(self):
        self.result = {}
        # self.result2 ={}
        self.__dp = 4
        self.list1=[]#frc31 > threshold
        self.list2=[]#frc32 > threshold
        self.info={}
        # self.info2={}

    def read_result(self,path,dp=4):
        import numpy as np
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



