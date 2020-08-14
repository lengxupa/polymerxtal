from abc import ABC, abstractmethod

class Observer(ABC):

    @abstractmethod
    def update(self, subject):
        pass


class XYZ_observer(Observer):

    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print('Printing xyz subject.xyz')


class energy_observer(Observer):

    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print('Printing energy', subject.energy)
