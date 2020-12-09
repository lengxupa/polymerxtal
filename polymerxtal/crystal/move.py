import math
import numpy as np


class Particle:
    '''Base class for a Particle. Particles
    can be translated, rotated, reflected, etc'''
    def move(self, transformer, **kwargs):
        method_name = 'move_{}'.format(self.__class__.__name__.lower())
        try:
            visit = getattr(transformer, method_name)
        except AttributeError:

            print("WARNING: The transformer {} does not contain the method {} for type {}".format(
                transformer.__class__.__name__, method_name, self.__class__.__name__))
            return
        return visit(self, **kwargs)


class Sphere(Particle):
    '''A Sphere is a type of particle, for instance,
    point particles with Lennard-Jones potential'''
    def __init__(self, center):
        self.center = center


class Gay_Berne(Particle):
    '''A Gay Berne particle is an anisotropic model
    of the 12-6 Lennard-Jones potential. It is used
    exstensively to model liquid crystals'''
    def __init__(self, center, vector):
        self.center = center
        self.vector = vector


class Cluster(Particle):
    '''A cluster is a collection of particles. It can
    be composed of spherical, anisotropic, patchy or
    any other particle'''
    def __init__(self):
        self.particles = []

    def add_particle(self, particle):
        self.particles.append(particle)


class Translator:
    def move_sphere(self, sphere, **kwargs):
        for key in kwargs:
            if key not in ['array']:
                raise KeyError('Translator - move_sphere: Must input array')

        array = kwargs['array']

        #print('Translating sphere should be very straightforward')
        sphere.center[0] += array[0]
        sphere.center[1] += array[1]
        sphere.center[2] += array[2]

    def move_gay_berne(self, gay_berne, **kwargs):
        for key in kwargs:
            if key not in ['array']:
                raise KeyError('Translator - move_gay_berne: Must input array')

        array = kwargs['array']

        #print('Translating solid should involve COM computation and subsequent translation')
        gay_berne.center[0] += array[0]
        gay_berne.center[1] += array[1]
        gay_berne.center[2] += array[2]

    def move_cluster(self, cluster, **kwargs):
        for key in kwargs:
            if key not in ['array']:
                raise KeyError('Translator - move_cluster: Must input array')

        array = kwargs['array']

        #print('Translating cluster should involve COM computation and subsequent translation')
        for particle in cluster.particles:
            particle.move(self, **kwargs)


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


class Rotator:
    def move_sphere(self, sphere, **kwargs):
        #print('Rotating sphere is nonsense due to rotation invariance')
        pass

    def move_gay_berne(self, gay_berne, **kwargs):
        for key in kwargs:
            if key not in ['axis', 'theta']:
                raise KeyError('Rotator - move_gay_berne: Must input axis and theta')

        axis = kwargs['axis']
        theta = kwargs['theta']

        matrix = rotation_matrix(axis, theta)

        #print('Rotating Gay Berne using Eulerian angles')
        #gay_berne.vector[0] = 1
        gay_berne.vector = np.dot(matrix, gay_berne.vector)

    def move_cluster(self, cluster, **kwargs):
        for key in kwargs:
            if key not in ['axis', 'theta']:
                raise KeyError('Rotator - move_cluster: Must input axis and theta')

        axis = kwargs['axis']
        theta = kwargs['theta']

        matrix = rotation_matrix(axis, theta)
        #print('Rotating cluster using Eulerian angles')
        for particle in cluster.particles:
            particle.move(self, **kwargs)
            particle.center = np.dot(matrix, particle.center)


def main():

    translator = Translator()
    rotator = Rotator()

    origin = [0.0, 0.0]
    vector = [1.0, 0.0]

    argon = Sphere(origin)
    krypton = Sphere(origin)
    liquid_crystal = Gay_Berne(origin, vector)

    cluster = Cluster()
    cluster.add_particle(argon)
    cluster.add_particle(krypton)
    cluster.add_particle(liquid_crystal)

    argon.move(translator)
    argon.move(rotator)

    krypton.move(translator)
    krypton.move(rotator)

    liquid_crystal.move(translator)
    liquid_crystal.move(rotator)

    cluster.move(translator)
    cluster.move(rotator)


if __name__ == "__main__":
    main()
