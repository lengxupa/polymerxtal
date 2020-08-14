from abc import abstractmethod, ABC
import numpy as np
import MDAnalysis
import mdtraj as md

class TrajectoryAdapter(ABC):
    @abstractmethod
    def compute_center_of_mass():
        pass
        
    @abstractmethod
    def compute_radius_of_gyration():
        pass


class MDAnalysisAdapter(TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = MDAnalysis.Universe(filename)
        print('Selected MDAnalysis.')
	
    def compute_center_of_mass(self):
        mass_by_frame = np.ndarray(shape=(len(self.trajectory.trajectory), 3))
        for ts in self.trajectory.trajectory:
            mass_by_frame[ts.frame] = self.trajectory.atoms.center_of_mass(compound='segments')
        return mass_by_frame
	
    def compute_radius_of_gyration(self):
        rg_by_frame = np.empty(len(self.trajectory.trajectory))
        for ts in self.trajectory.trajectory:
            rg_by_frame[ts.frame] = self.trajectory.atoms.radius_of_gyration()
        return rg_by_frame
		
		
class MDTrajAdapter(TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = md.load_pdb(filename)
        print('Selected MDTraj.')
	
    def compute_center_of_mass(self):
        return 10 * md.compute_center_of_mass(self.trajectory)
	
    def compute_radius_of_gyration(self):
        return 10 * md.compute_rg(self.trajectory)