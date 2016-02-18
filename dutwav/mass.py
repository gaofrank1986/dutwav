class Mass(object):
    """
      class for store mass information
    """
    def __init__(self):
        self._flag = -1
        self.rot_ctr = None
        self.mass_ctr = None
        self.water_plane_info = None
        self.displacement = None
        self.disp_ctr = None
        self.mass = None
        self.mass_mtx = None
        self.hyres_mtx = None
        self.lin_damp_mtx = None
        self.quad_damp_mtx = None
        self.stif_mtx = None
        self.info = {}

    def read_mass_info(self,path):
        import numpy as np
        with open(path,"r") as f:
            tmp = f.readline().strip()
            flag =[int(i) for i in f.readline().split()][0]
            self.info[tmp] = flag
            self._flag = flag

            tmp = f.readline().strip()
            rotation_ctr =[float(i) for i in f.readline().split()]
            self.info[tmp] = rotation_ctr
            self.rot_ctr = np.array(rotation_ctr)

            tmp = f.readline().strip()
            mass_ctr =[float(i) for i in f.readline().split()]
            self.info[tmp] = mass_ctr
            self.mass_ctr = np.array(mass_ctr)

            tmp = f.readline().strip()
            water_plane_info =[float(i) for i in f.readline().split()]
            self.info[tmp] = water_plane_info
            self.water_plane_info = water_plane_info

            tmp = f.readline().strip()
            displace =[float(i) for i in f.readline().split()]
            self.info['displacement/volume'] = displace[0]
            self.info['displacement/buoy center'] = displace[1:4]
            self.displacement = displace[0]
            self.disp_ctr = np.array(displace[1:4])

            tmp = f.readline().strip()
            mass =[float(i) for i in f.readline().split()][0]
            self.info['mass'] = mass
            self.mass = mass
            
            mat_mass = []
            for i in range(3):
                mat_mass.append([float(i) for i in f.readline().split()])
            mat_mass = np.array(mat_mass)
            self.info['mass matrix'] = mat_mass
            self.mass_mtx = mat_mass
            mat_list = []
            for j in range(4):
                tmp = f.readline()
                mat = []
                for i in range(6):
                    mat.append([float(i) for i in f.readline().split()])
                mat = np.array(mat)
                mat_list.append(mat)
            self.info['hydro_storing matrix'] = mat_list[0]
            self.info['stiffness matrix'] = mat_list[1]
            self.info['linear viscous damping matrix'] = mat_list[2]
            self.info['quadratic viscous damping matrix'] = mat_list[3]
            self.hyres_mtx = mat_list[0]
            self.stif_mtx = mat_list[1]
            self.lin_damp_mtx = mat_list[2]
            self.quad_damp_mtx = mat_list[3]
