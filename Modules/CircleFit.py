import numpy as np
import matplotlib.pyplot as plt
import Modules.LinearFit as lf


class CircleFit:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y
        self.Q = 0
        x_m = np.mean(x)
        y_m = np.mean(y)

        # calculation of the reduced coordinates
        u = (x) - x_m
        v = (y) - y_m

        # linear system defining the center (uc, vc) in reduced coordinates:
        #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
        #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
        Suv  = np.sum(u*v)
        Suu  = np.sum(u**2)
        Svv  = np.sum(v**2)
        Suuv = np.sum(u**2 * v)
        Suvv = np.sum(u * v**2)
        Suuu = np.sum(u**3)
        Svvv = np.sum(v**3)

        # Solving the linear system
        A = np.array([ [ Suu, Suv ], [Suv, Svv]])
        B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
        uc, vc = np.linalg.solve(A, B)

        self.xc_1 = x_m + uc
        self.yc_1 = y_m + vc

        # Calculates distances to centre (xc_1, yc_1)
        self.Ri_1     = np.sqrt(((x)-self.xc_1)**2 + ((y)-self.yc_1)**2)
        self.R_1      = np.mean(self.Ri_1)
        self.residu_1 = np.sum((self.Ri_1-self.R_1)**2)

        self.p = self.R_1 *300.*4 ###### B=20T !!!!!!!!!!!!!!!!!!!!!!!!!!!


        # calculate angle of tangent to x axis, phi
        m_radius = self.yc_1/self.xc_1
        m_tangent = -1/m_radius

        # get the charge in order to manipulate phi value
        particle_dir = np.array([x[0], y[0], 0])
        circle_dir = np.array([self.xc_1, self.yc_1, 0])
        cross = np.cross(particle_dir, circle_dir)
        if cross[2]>0:
            self.Q = -1
        elif cross[2]<0:
            self.Q = 1
        else:
            self.Q = 0

        self.phi = np.arctan2(self.yc_1, self.xc_1) + (np.pi/2)*self.Q
        while self.phi > np.pi: self.phi -= 2*np.pi
        while self.phi < -np.pi: self.phi += 2*np.pi

    @property
    def Q(self):
        return self._Q
    @Q.setter
    def Q(self, value):
        self._Q = value

        
    def __str__(self):
        output = f"Center (x,y) = ({self.xc_1}, {self.yc_1}), Radius = {self.R_1}, Residuals = {self.residu_1}"
        return output

    def plot(self, xlabel = 'x / m', ylabel = 'y / m', color = 'm', label = 'Data', marker = 'o', filename = 'unnamedcirclefit.png'):

        # tangent parameters
        m_radius = self.yc_1/self.xc_1
        m_tangent = -1/m_radius
        c_tangent = 0
        x_tangent = np.linspace(0,self.x[np.argmax(abs(self.x))], 100)

        # plot
        plt.figure(figsize = (8, 8))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.errorbar(self.x, self.y,
                     color = color, marker = marker, linestyle = '', label = label)
        circle = plt.Circle((self.xc_1, self.yc_1), radius=self.R_1, color='black', fill=False)
        plt.gca().add_artist(circle)
        plt.plot(x_tangent, m_tangent*x_tangent+c_tangent, linestyle='--', color='c', label='Tangent')
        plt.axhline(y=0, linestyle='-', color='g')
        #plt.xlim([ 0, self.x[np.argmax(abs(self.x))] ])
        plt.legend()
        plt.savefig(filename)

        return plt

    def get_radius(self):
        return self.R_1

    def get_tot_momentum(self):
        return self.p

    def get_angle_phi(self):
        return self.phi

    def get_px(self):
        return self.p*np.cos(self.phi)

    def get_py(self):
        return self.p*np.sin(self.phi)

    def get_Q(self):
        return self.Q

