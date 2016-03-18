from numpy import sqrt
class Wave(object):
    import numpy as np
    k = np.nan
    freq = np.nan
    wav_len = None
    A = np.nan
    beta = None
    d = None 

    def compute_wav_num(sigma,h):
        if (sigma < 0):
            print "Error,in compute_wav_num, wave_frequency = ", sigma
        else:
            if (h > 0):
                B = _g*h
                Y = sigma**2*h/_g
                A = 1.0/(1.0+Y*(0.66667+Y*(0.35550+Y*(0.16084+Y*(0.063201+Y*(0.02174+Y*(0.00654+Y*(0.00171+Y*(0.00039+Y*0.00011)))))))))
                _wave_celerity=SQRT(B/(Y+A))
                _wave_num = sigma/_wave_celerity
            else:
                if(h < 0):
                    _wav_num = sigma**2/_g
        return _wav_num

    def read_data(self,path):
        G = 9.807

        with open(path,"r") as f:
            data1 = [int(i) for i in f.readline().split()]
            data2 = [float(i) for i in f.readline().split()]
        
        self.d = data2[0]
        self.A = data2[1]
        self.beta = data2[3]

        if (data1[1] ==0):
            self.k = data2[2]
            if self.d < 0:
                self.freq = sqrt(G*self.k)
            else:
                self.freq = sqrt(G*self.k*tanh(self.k*self.d))
        else:
            self.wave_len = data2[2]
            self.k = compute_wav_num(self.freq,self.d)

        # self.A = data2[1]
    @property
    def info(self):
        print "wave number: ",self.k,"\n"
        print "wave frequency: ",self.freq,"\n"
        print "wave amplitude: ",self.A,"\n"
        print "water depth: ",self.d,"\n"
        print "wave angle:",self.beta,"\n"

