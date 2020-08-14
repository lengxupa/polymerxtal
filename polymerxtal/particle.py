from abc import ABC, abstractmethod

class Particle:

    '''Base class for a Particle. Particles
    can be translated, rotated, reflected, etc'''

    def move(self, transformer):
        method_name = 'move_{}'.format(self.__class__.__name__.lower())
        try:
            visit = getattr(transformer, method_name)
        except AttributeError:

            print("WARNING: The transformer {} does not contain the method {} for type {}".format(transformer.__class__.__name__, method_name, self.__class__.__name__))
            return
        return visit(self)

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
