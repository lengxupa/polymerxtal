
from abc import ABC, abstractmethod
import numpy as np


class Potential(ABC):
   @abstractmethod
   def get_energy(self):
       pass


class LJ(Potential):
   def __init__(self, **kwargs):
       for key in kwargs:
           if key not in ['epsilon', 'sigma']:
               raise KeyError('LJ potential: Must input epsilon and sigma')

       self.sigma = kwargs['sigma']
       self.epsilon = kwargs['epsilon']

   def get_energy(self, r):
       return 4 * self.epsilon * ((self.sigma / r)**12 - (self.sigma / r)**6)


class Buckingham(Potential):
   def __init__(self, **kwargs):
       for key in kwargs:
           if key not in ['A', 'C', 'rho']:
               raise KeyError('Buckingham potential: Must input A, C and rho')

       self.rho = kwargs['rho']
       self.A = kwargs['A']
       self.C = kwargs['C']

   def get_energy(self, r):
       return self.A * np.exp(-r / self.rho) - self.C / r**6

