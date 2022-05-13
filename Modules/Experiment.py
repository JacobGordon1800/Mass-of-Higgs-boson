import numpy as np

class Experiment:
    '''
    Class to define the Geometry of an experiment
    initialisation may take two optional arguments:
    * featuresize that is the smallest element size (in m)
    * rng is a numpy random number generator
       
    An arbitrary number of experimental features are added via add_volume like the class Cuboid coded above
    
    The class has four non-trivial methods:
    * scan_volume_change: scan the straight line between two points for a change of volume
    * do_eloss: take away energy from a particle according to material properties
    * do_mult_scatter: change the direction of a particle according to multiple scattering
    * detect_particles: simple particle position 'detection'
    '''

    def __init__(self, featuresize = 1E-4, rng = None):
        self.volumes = [] # an empty list of volumes
        self.minfeaturesize = featuresize
        
        # possibility to use the same random number generator everywhere
        # if argument is skipped, this creates its own instance of the default random number generator
        if rng is None:
            self._rng = np.random.default_rng()
        else:
            self._rng = rng

    def add_volume(self, shape):
        print(f"Adding shape")
        print(shape)
        self.volumes.append(shape)
        return len(self.volumes)

    def get_volume(self, xyz):
        for id in range(len(self.volumes)):
            if self.volumes[id].is_in_volume(xyz):
                return id

        # if we arrived here, we are outside everything
        return -1

    def scan_volume_change(self, oldpos, delta, lastVolume):
        # scan the straight line between last (oldpos) along step (delta)
        # for a missed feature of the experiment

        dist = np.sqrt(np.sum(delta**2)) # this is the distance in x-y-z
        nsteps = int(dist/self.minfeaturesize) + 1
        if nsteps <= 1:
            # step was small enough, continue
            return 1.

        for n in range(1,nsteps):
            if self.get_volume(oldpos+(n/nsteps)*delta) != lastVolume:
                if n == 1:
                    return (0.5 / nsteps)
                else:
                    return (n-1.0) / nsteps
        return 1.

    def do_eloss(self, particle, dist):
        '''
        Method checks if particle is in a material and if so, reduces its energy by dE/dx*distance
        '''
        volume = self.get_volume(particle.xyz)
        if volume >= 0 and self.volumes[volume].material != 0:
            lostE = self.volumes[volume].material.get_energy_loss(particle)*dist
            particle.reduce_energy(lostE)

    def do_mult_scatter(self, particle, dist):
        '''
        Method checks if particle is in a material and if so, scatters it by theta0
        '''
        volume = self.get_volume(particle.xyz)
        if volume >= 0 and self.volumes[volume].material != 0:
            theta0 = self.volumes[volume].material.get_theta0(particle, dist)
            particle.apply_small_rotation(self._rng.standard_normal(2)*theta0)

    def detect_particles(self, track):
        '''
        Method takes a track and checks crossings with all volumes
        '''

        detection = np.zeros((len(self.volumes), 8))
        ncross = np.zeros(len(self.volumes))

        # loop over all points and average over those matching volumes
        for ipoint in range(len(track)):
            idet = self.get_volume(track[ipoint][5:8])
            if idet >= 0:
                detection[idet] += track[ipoint]
                ncross[idet] += 1

        for idet in range(len(self.volumes)):
            if ncross[idet] > 0:
                detection[idet] /= ncross[idet]

        return detection
    
