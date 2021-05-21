from abc import ABC, abstractmethod
import numpy as np


class Subject:
    def __init__(self):
        self._observers = []

    def attach(self, observer):
        if observer not in self._observers:
            self._observers.append(observer)

    def detach(self, observer):
        if observer in self._observers:
            self._observers.remove(observer)

    def notify(self):
        for observer in self._observers:
            observer.update(self)


class Observer(ABC):
    @abstractmethod
    def update(self, subject):
        pass


class System(Subject):
    def __init__(self, mol_number=0):
        Subject.__init__(self)
        self._xyz = np.zeros((mol_number, 3))
        self.current_step = 0
        self.compute_energy()

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        self._xyz = xyz
        self.current_step += 1
        self.compute_energy()
        self.notify()

    def compute_energy(self):
        # Ideally we would implement a full energy computation
        self.energy = np.random.rand()


class XYZ_observer(Observer):
    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print("Printing xyz subject.xyz")


class Energy_observer(Observer):
    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print("Printing energy", subject.energy)


def main():
    # Create the system with a desired number of particles
    mc_system = System(mol_number=15)

    xyz_observer = XYZ_observer(print_freq=2)
    energy_observer = Energy_observer()

    mc_system.attach(xyz_observer)
    mc_system.attach(energy_observer)

    total_mcsteps = 5
    for this_mcstep in range(total_mcsteps):

        # Ideally you would perform some MC moves here
        # that alter the state of the system and
        # don't violate detailed balance!
        mc_system.xyz += 0.1


if __name__ == "__main__":
    main()
