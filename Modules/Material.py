import Modules.Particle as p
import numpy as np

class Material:
    '''
    A class to collect material-related effects
    such as Ionisation energy loss and Multiple scattering
    '''

    # some constants
    K = 0.0307075 # MeV m^2 kg^-1 (Eloss constant)
    me = 0.511   # MeV (electron mass)

    def __init__(self, rho, Z, A):
        self.rho = rho
        self.Z = Z
        self.A = A

        self.I = 0.0000135 * Z

        if rho <= 0. or Z <= 0. or A <= 0.:
            self.X0 = -1.
        else:
            self.X0 = 7164*A/(rho*Z*(Z+1)*np.log(287/np.sqrt(Z)))

    def __str__(self):
        output = f"Material with density {self.rho} kg/m^3, Z = {self.Z}, A = {self.A}, X0 = {self.X0}m"
        return output

    def get_energy_loss(self, particle):
        '''
        Ionisation energy loss in MeV/m
        '''
        if self.X0 < 0.:
            return 0.

        b = particle.beta()
        g = particle.gamma()
        z = particle.Q
        M = particle.m;
        meM = self.me/M;
        Wmax = 2*self.me*b*b*g*g/(1+2*g*meM+meM*meM)

        # textbook formula in MeV/cm
        eLoss = self.K*z*z*self.rho*self.Z/self.A/(b*b)*(0.5*np.log(2*self.me*b*b*g*g*Wmax/(self.I*self.I)) - b*b)
        return eLoss

    def get_theta0(self, particle, depth):
        '''
        Scattering angle theta0 when traversing material of depth
        '''

        if self.X0 < 0.:
            return 0.

        b = particle.beta()
        z = particle.Q
        p = particle.momentum()
        return 13.6/(b*p)*abs(z)*np.sqrt(depth/self.X0)*(1+0.038*np.log(depth/self.X0))
    
