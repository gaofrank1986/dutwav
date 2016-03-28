
from dutwav.mesh import Mesh
from numpy import sign

class Box(Mesh):
    def __init__(self):
        Mesh.__init__(self)
        self.x1=-0.5
        self.x2=0.5
        self.y1=0.5
        self.y2=-0.5
        self.z=-1.

    def box_body_renrml(self):
        
        self._mark_elems_at(0,self.x1,'x1')
        self._mark_elems_at(0,self.x2,'x2')
        self._mark_elems_at(1,self.y1,'y1')
        self._mark_elems_at(1,self.y2,'y2')
        self._mark_elems_at(2,self.z,'z')
        mx1=self.extract_mesh(['x1'])
        mx2=self.extract_mesh(['x2'])
        my2=self.extract_mesh(['y2'])
        my1=self.extract_mesh(['y1'])
        mz=self.extract_mesh(['z'])
        mx1._update_all_nrml([sign(self.x1),0.,0.0])
        mx2._update_all_nrml([sign(self.x2),0.,0.0])
        my1._update_all_nrml([0.,sign(self.y1),0.0])
        my2._update_all_nrml([0.,sign(self.y2),0.0])
        mz._update_all_nrml([0.,0.0,sign(self.z)])
        mx1.devour_mesh(mx2)
        mx1.devour_mesh(my2)
        mx1.devour_mesh(my1)
        mx1.devour_mesh(mz)
        return mx1


class Cylinder(Mesh):
    def __init__(self):
        Mesh.__init__(self)
        self.z1=-1
        self.z2=-6
        self.r1=1
        self.r2=8

    def cylinder_nrml(self):

        self._mark_elems_at(2,self.z1,'z1')
        self._mark_elems_at(2,self.z2,'z2')
        self._mark_elem_withr(self.r1-0.1,self.r1+0.1,'r1')
        self._mark_elem_withr(self.r2-0.5,self.r2+0.5,'r2')
        self._count_elem()
        mz1=self.extract_mesh(['z1'])
        mz2=self.extract_mesh(['z2'])
        mr2=self.extract_mesh(['r2'])
        # mr2.draw_model()
        mr1=self.extract_mesh(['r1'])
        mr1._update_tag(['extract'],'body')
        mr2._update_tag(['extract'],'body')
        mr1._redo_bd_nrml()
        mr1._reverse_nrml()
        mr2._redo_bd_nrml()
        mz1._update_all_nrml([0.,0.0,-sign(self.z1)])
        mz2._update_all_nrml([0.,0.0,sign(self.z1)])
        mr2.devour_mesh(mz2)
        mr2.devour_mesh(mr1)
        mr2.devour_mesh(mz1)
        return mr2






