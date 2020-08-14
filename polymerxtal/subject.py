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


class System(Subject):

    def __init__(self, mol_number=0):
        Subject.__init__(self)

    # Initialize particle positions
        self._xyz = np.zeros((mol_number, 3))

        self.current_step = 0

    # Compute initial energy
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


