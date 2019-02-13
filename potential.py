from abc import ABC, abstractmethod
import numpy as np


class Potential(ABC):
    @abstractmethod
    def get_energy():
        pass


# Needs docstrings
class LJ(Potential):
    def __init__(self, kwargs):

        for key, value in kwargs.items():
            if key not in ["epsilon", "sigma"]:
                raise KeyError("LJ potential: Must input epsilon and sigma")

        self.sigma = kwargs["sigma"]
        self.epsilon = kwargs["epsilon"]

    def get_energy(self, r):
        return 4 * self.epsilon * ((self.sigma / r)**12 - (self.sigma / r)**6)


# Needs docstrings
class Buckingham(Potential):
    def __init__(self, kwargs):

        for key, value in kwargs.items():
            if key not in ["A", "C", "rho"]:
                raise KeyError("Buckingham potential: Must input A, C and rho")

        self.rho = kwargs["rho"]
        self.A = kwargs["A"]
        self.C = kwargs["C"]

    def get_energy(self, r):
        return self.A * np.exp(-r / self.rho) - self.C / r**6


def potential_factory_v1(potential_type, **kwargs):

    if potential_type == "LJ":
        return LJ(kwargs)

    elif potential_type == "Buckingham":
        return Buckingham(kwargs)

    else:
        raise Exception("Potential type not found")


def potential_factory_v2(potential_type, **kwargs):
    cls_dict = dict(LJ=LJ, Buckingham=Buckingham)

    if potential_type not in cls_dict.keys():
        raise Exception("Potential type not found")

    cls = cls_dict[potential_type]

    return cls(kwargs)


pt = potential_factory_v1("Buckingham", A=4.0, rho=10.0, C=10)
print("Buckingham     v1: ", pt.get_energy(r=10.0))

pt = potential_factory_v1("LJ", epsilon=40.0, sigma=3.54)
print("Lennard-Jones  v1: ", pt.get_energy(r=10.0))

pt = potential_factory_v2("Buckingham", A=4.0, rho=10.0, C=10)
print("Buckingham     v2: ", pt.get_energy(r=10.0))

pt = potential_factory_v2("LJ", epsilon=40.0, sigma=3.54)
print("Lennard-Jones  v2: ", pt.get_energy(r=10.0))
