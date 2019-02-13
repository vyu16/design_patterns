import mdtraj as mt
import MDAnalysis as ma
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
        print("Using MDTraj")

    def compute_center_of_mass(self):
        return 10 * mt.compute_center_of_mass(self.trajectory)

    def compute_radius_of_gyration(self):
        return 10 * mt.compute_rg(self.trajectory)


class MDAnalysisAdapter(TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = ma.Universe(filename)
        print("Using MDAnalysis")



mda = MDTrajAdapter("protein.pdb")

print(mda.compute_center_of_mass())
print(mda.compute_radius_of_gyration())
