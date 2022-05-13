#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np

# units are MeV (with supressed factors of c) and metres, seconds
# charge is in units of elementary charge

class Particle:
    c = 299792458. # m/s

    def __init__(self, E = -1., px=0., py=0., pz=0.,
                 t = 0., x = 0., y = 0., z = 0., m = -1, Q = -1):
        self._Epxpypztxyz = np.zeros(8)

        self.px = px
        self.py = py
        self.pz = pz
        self.t = t
        self.x = x
        self.y = y
        self.z = z

        self.Q = Q

        # optional set either E or mass
        if E < 0.:
            self.E = np.sqrt(self.momentum()**2 + m*m)
            self.m = m
        else:
            self.E = E
            self.m = np.sqrt(E*E - self.momentum()**2)

    @property
    def E(self):
        return self._Epxpypztxyz[0]
    @E.setter
    def E(self, value):
        self._Epxpypztxyz[0] = value
    @property
    def px(self):
        return self._Epxpypztxyz[1]
    @px.setter
    def px(self, value):
        self._Epxpypztxyz[1] = value
    @property
    def py(self):
        return self._Epxpypztxyz[2]
    @py.setter
    def py(self, value):
        self._Epxpypztxyz[2] = value
    @property
    def pz(self):
        return self._Epxpypztxyz[3]
    @pz.setter
    def pz(self, value):
        self._Epxpypztxyz[3] = value

    @property
    def t(self):
        return self._Epxpypztxyz[4]
    @t.setter
    def t(self, value):
        self._Epxpypztxyz[4] = value
    @property
    def x(self):
        return self._Epxpypztxyz[5]
    @x.setter
    def x(self, value):
        self._Epxpypztxyz[5] = value
    @property
    def y(self):
        return self._Epxpypztxyz[6]
    @y.setter
    def y(self, value):
        self._Epxpypztxyz[6] = value
    @property
    def z(self):
        return self._Epxpypztxyz[7]
    @z.setter
    def z(self, value):
        self._Epxpypztxyz[7] = value

    @property
    def m(self):
        return self._m
    @m.setter
    def m(self, value):
        self._m = value
        
    @property
    def Q(self):
        return self._Q
    @Q.setter
    def Q(self, value):
        self._Q = value
    @property
    def state(self):
        return (self._Epxpypztxyz)
    @property
    def xyz(self):
        return (self._Epxpypztxyz[5:8])
    #@xyz.setter
    #def xyz(self, value):
    #    self._Epxpypztxyz[5:8] = value
        
    def __str__(self):
        output = f"Particle with mass {self.m} MeV and charge {self.Q:+}e\n"
        output += f"(t,x,y,z) = ({self.t}s, {self.x}m, {self.y}m, {self.z}m)\n"
        output += f"(E,px,py,pz) = ({self.E}, {self.px}, {self.py}, {self.pz}) MeV"
        return output

    def momentum(self):
        return np.sqrt(self.px**2 + self.py**2 + self.pz**2)

    def mass(self):
        return self.m
        
    def gamma(self):
        return (self.E/self.m)

    def beta(self):
        return (self.momentum()/self.E)

    def radius(self, B):
        return (self.momentum()/(self.c*1e-6*B))

    def __add__(self, other):
        sum = self._Epxpypztxyz + other._Epxpypztxyz
        return Particle(sum[0], sum[1], sum[2], sum[3], sum[4], sum[5], sum[6], sum[7],
                        0., self.Q+other.Q)

    def apply_small_rotation(self, dtheta_xz_yz):
        '''
        calculate new momentum direction
        approximates this into two sequential x-z and y-z changes
        '''
        pxz = np.sqrt(self.px*self.px + self.pz*self.pz);
        theta = np.arctan2(self.px, self.pz) + dtheta_xz_yz[0]
        self.px = pxz*np.sin(theta)
        self.pz = pxz*np.cos(theta)

        pyz = np.sqrt(self.py*self.py + self.pz*self.pz);
        theta = np.arctan2(self.py, self.pz) + dtheta_xz_yz[1]
        self.py = pyz*np.sin(theta)
        self.pz = pyz*np.cos(theta)

    def reduce_energy(self, Eloss):
        '''reduce energy of particle while keeping direction the same'''

        Enew = self.E - Eloss
        if Enew <= self.m:
            # all kinetic energy lost, put particle at rest
            self.E = self.m
            self._Epxpypztxyz[1:4] = 0.
        else:
            pnew = np.sqrt(Enew*Enew - self.m*self.m)
            factor = pnew/self.momentum()
            self.E = Enew
            self._Epxpypztxyz[1:4] *= factor

