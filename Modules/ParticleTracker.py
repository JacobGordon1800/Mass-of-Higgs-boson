import numpy as np
import copy
import Modules.Particle as p

class ParticleTracker:
    '''
    Particle Tracker class
    '''

    def __init__(self, particle, tfinal, dtmax):
        # deepcopy ensures to really create a new object with the same content
        # a simple "=" opertion would just create a reference to the same object
        self.particle_in = copy.deepcopy(particle)
        self.particle_out = copy.deepcopy(particle)
        
        # the maximum step size and the target time
        self._dtmax = dtmax
        self._tfinal = tfinal

        # step counter
        self._nsteps = 0
        
        # prepare storage for the E-momentum-time-pos (8 components) as the particle moves along
        # calculate the expected number of steps and generously double
        nsteps = self._tfinal/self._dtmax
        self.track = np.zeros( (int(nsteps*2)+1,8) )
        # save the first point
        self.track[0] = self.particle_in.state[:]

        # E + B field at "current" location of the particle, start at all zero
        self._Efield = np.zeros(3)
        self._Bfield = np.zeros(3)

    def propagate(self, experiment):

        # initialise first step
        t = self.particle_in.t
        
        while t < self._tfinal:

            # take a step by dt, make the last step reach the target
            dt = self._dtmax
            if t+1.05*dt > self._tfinal:
                dt = self._tfinal - t

            lastVolume = experiment.get_volume(self.particle_out.xyz)
            lastStep = copy.deepcopy(self.particle_out)
            
            # !large debug output!, this should be commented in "production runs"
            #print(f"Step {self._nsteps} with size {dt} from volume {lastVolume}")
            #print(self.particle_out)
            # end of debug output

            step = self.dostep(dt)

            # check, if the last step crossed into or over a different experimental volume
            reduceStep = experiment.scan_volume_change(lastStep.xyz, step[5:8], lastVolume)
            if reduceStep < 1.0:
                # yes, we crossed a boundary: go back and move forward with a reduced step size
                self.particle_out = lastStep
                step = self.dostep(dt*reduceStep)
            
            # implement Multiple Coulomb scattering
            dist_xyz = np.sqrt(np.sum((step[5:8])**2)) # this is the distance in x-y-z
            experiment.do_mult_scatter(self.particle_out, dist_xyz)

            # implement Ionisation Energy Loss
            experiment.do_eloss(self.particle_out, dist_xyz)

            # Abort if particle stopped moving
            if self.particle_out.momentum() <= 0.:
                print("Particle out of energy, done.")
                break

            # !large debug output!, this should be commented in "production runs"
            #print(f"After step {self._nsteps} in volume {lastVolume} dist {dist_xyz}")
            #print(self.particle_out)
            # end of debug output


            # update time and number of steps taken
            t = self.particle_out.t
            self._nsteps += 1
            
            # store the current E-momentum-time-pos
            # make sure there's enough space in the track array
            if self.track.shape[0] <= self._nsteps:
                extension = np.zeros((self._nsteps, 8))
                self.track = np.concatenate((self.track,extension), axis=0)
            self.track[self._nsteps] = self.particle_out.state[:]
            

        # we've reached the target time
        # cut the track data and return the number of taken steps
        self.track = self.track[:self._nsteps+1]
        return self._nsteps


    def update_EBfields(self, pos):
        '''
        define vectorial E (V/m) and B (Tesla) fields at position 
        given by pos (that has pos[4]=t, pos[5]=x, pos[6]=y, pos[7]=z'''

        # extract the coordinates you need for nicer equations below
        x = pos[5]
        y = pos[6]
        z = pos[7]
        
        # all zeros unless changed, e.g. a constant field of 4T in z-direction for all space
        self._Bfield[2] = 4. ######################################## B=20T !!!!!!!!!!!!!!!!!!!!!!!

        return self._Efield, self._Bfield

    def dostep(self, dt):
        '''
        Make one Step of dt - Runge-Kutta 2nd order Method (could be replaced by Euler or Runge-Kutta4 Method)
        '''
        
        # get the current state of particle as array of 8 elements: y=(E,px,py,pz,t,x,y,z)
        y = self.particle_out.state

        # the two RK2 steps
        k1 = self.get_differential(y)
        y2 = y + dt*k1
        k2 = self.get_differential(y2)
        step = dt*0.5*(k1 + k2)

        y += step
        return step

    
    def get_differential(self, y):
        '''
        calculate the differential f(t, y(t))
        input is current state y = (E, px, py, pz, t, x, y, z)
        output is the time differential if this vector dy/dt = d/dt (E, px, py, pz, t, x, y, z)
        that is equal to (dE/dt, dpx/dt, dpy/dt, dpz/dt, 1, vx, vy, vz, )
        component 0 is the change in Energy q*(v . E)
        components 1-3 are the vectorial Lorentz force q*c(E + (v x B))
        component 4 is trivially dt/dt=1
        components 5-8 are the relativistic speed given by v=c*p/E
        unfortuantely some numpy vector operations are very slow, so cannot use np.cross(v,B)
        '''

        c = p.Particle.c
        Q = self.particle_out.Q*1E-6 # conversion of units requires 1E-6 factor
        v = (c/y[0])*y[1:4] # relativistic speed c/E*p

        # update the E and B fields at coordinates y=(E, px, py, pz, t, x, y, z)
        E, B = self.update_EBfields(y)
        
        f = np.array([Q*(np.dot(E,v)),
                      Q*c*(E[0] + v[1]*B[2]-v[2]*B[1]),
                      Q*c*(E[1] + v[2]*B[0]-v[0]*B[2]),
                      Q*c*(E[2] + v[0]*B[1]-v[1]*B[0]),
                      1.,
                      v[0],
                      v[1],
                      v[2]])

        return f

