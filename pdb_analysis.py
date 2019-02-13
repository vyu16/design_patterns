import mdtraj as mt
import MDAnalysis as ma
import numpy as np
from abc import abstractmethod, ABC


class TrajectoryAdapter(ABC):
    @abstractmethod
    def compute_center_of_mass():
        pass

    @abstractmethod
    def compute_radius_of_gyration():
        pass


class MDTrajAdapter(TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = mt.load_pdb(filename)
        print("\nUsing MDTraj")

    def compute_center_of_mass(self):
        return 10 * mt.compute_center_of_mass(self.trajectory)

    def compute_radius_of_gyration(self):
        return 10 * mt.compute_rg(self.trajectory)


class MDAnalysisAdapter(TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = ma.Universe(filename)
        print("\nUsing MDAnalysis")

    def compute_center_of_mass(self):
        mass_by_frame = np.ndarray(shape=(len(self.trajectory.trajectory), 3))

        for ts in self.trajectory.trajectory:
            mass_by_frame[ts.frame] = self.trajectory.atoms.center_of_mass(
                compound="segments")

        return mass_by_frame

    def compute_radius_of_gyration(self):
        rg_by_frame = np.empty(len(self.trajectory.trajectory))

        for ts in self.trajectory.trajectory:
            rg_by_frame[ts.frame] = self.trajectory.atoms.radius_of_gyration()

        return rg_by_frame


mda = MDTrajAdapter("protein.pdb")

print(mda.compute_center_of_mass())
print(mda.compute_radius_of_gyration())

mdb = MDAnalysisAdapter("protein.pdb")

print(mdb.compute_center_of_mass())
print(mdb.compute_radius_of_gyration())
