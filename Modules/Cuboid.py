import numpy as np

class Cuboid:
    def __init__(self, lowx, highx, lowy, highy, lowz, highz, material=0):
        self._lowx = lowx
        self._highx = highx
        self._lowy = lowy
        self._highy = highy
        self._lowz = lowz
        self._highz = highz
        self.material = material
        
    def is_in_volume(self, xyz):
        '''
        test if point xyz=(x,y,z) is inside the cuboid
        '''

        if xyz[0] < self._lowx: return False
        if self._highx < xyz[0]: return False
        if xyz[1] < self._lowy: return False
        if self._highy < xyz[1]: return False
        if xyz[2] < self._lowz: return False
        if self._highz < xyz[2]: return False

        return True
        
    def is_in_volume(self, xyz):
        '''
        test if point xyz=(x,y,z) is inside the cuboid
        '''

        if self._lowx < xyz[0] < self._highx         and self._lowy < xyz[1] < self._highy         and self._lowz < xyz[2] < self._highz :
            return True
        else:
            return False
    
    def __str__(self):
        output = f"Cuboid with edge low ({self._lowx}, {self._lowy}, {self._lowz})"
        output += f" and edge high ({self._highx}, {self._highy}, {self._highz})\n"
        output += self.material.__str__()
        return(output)
