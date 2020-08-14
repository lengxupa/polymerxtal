import numpy as np
import MDAnalysis
import mdtraj as md

class MDAnalysisTrajectoryHandler:
    def __init__(self, filename):
        self.trajectory = MDAnalysis.Universe(filename)
	
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
		
		
class MDTrajTrajectoryHandler:
    def __init__(self, filename):
        self.trajectory = md.load_pdb(filename)
	
    def compute_center_of_mass(self):
        return 10 * md.compute_center_of_mass(self.trajectory)
	
    def compute_radius_of_gyration(self):
        return 10 * md.compute_rg(self.trajectory)


class TrajectoryAnalyzer:
    def __init__(self, filename):
        self.mda = MDAnalysisTrajectoryHandler(filename)
        self.mdt = MDTrajTrajectoryHandler(filename)

    def analysis(self):
        results = {}
        results['MDAnalysis_center_of_mass'] = self.mda.compute_center_of_mass()
        results['MDTraj_radius_of_gyration'] = self.mdt.compute_radius_of_gyration()
        return results

