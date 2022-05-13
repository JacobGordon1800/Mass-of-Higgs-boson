#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np
import Modules.Material as m
class Cylinder:
    def __init__(self, lowr, highr, lowz, highz, material = 0):
        self._lowr = lowr
        self._highr = highr
        self._lowz = lowz
        self._highz = highz
        self.material = material

    def is_in_volume(self, xyz):
        r = np.sqrt(xyz[0]**2 + xyz[1]**2)
        if r > self._highr: return False
        if r < self._lowr: return False
        if xyz[2] < self._lowz: return False
        if xyz[2] > self._highz: return False
        return True

    def __str__(self):
        output = f"Cylinder with radius ({self._lowr} to {self._highr})"
        output += f"with high height ({self._highz})\n"
        output += f" and low height ({self._lowz})\n"
        output += self.material.__str__()
        return output


# In[11]:

'''
#setting up cylinder

x = 0.5
y = 0.5
r = np.sqrt(x**2 + y**2)
lowz = -0.5
highz = 0.5
ironcylinder = Cylinder(r = r, lowz = lowz, highz = highz,
               material = m.Material(rho = 7870, Z = 26, A = 55.845))
print(ironcylinder)

#testing 
rtest = np.linspace(r-.01, r+.01, 3)
ztest = np.linspace(lowz-.01, highz+.01, 3)

for r in rtest:
    for z in ztest:
        print(f"({r:.5}, {z:.5}) is within the iron cylinder: ", ironcylinder.is_in_volume((r,z)))
        
'''
