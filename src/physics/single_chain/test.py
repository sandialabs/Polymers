"""Module to test the local module.

"""
import numpy as np


def integrate(function, lower_lim, upper_lim, num_points):
    """Function for integration.

    """
    d_x = (upper_lim - lower_lim)/num_points
    return function(
        lower_lim + (0.5 + np.arange(0, num_points))*d_x
    ).sum()*d_x


class BasicParameters:
    """Class for basic testing parameters.

    """
    def __init__(self):
        self.abs_tol = 1e-7
        self.rel_tol = 1e-5
        self.log_log_tol = 1e-1
        self.log_log_scale = 12e-1
        self.number_of_loops = 8
        self.zero = 1e-6
        self.points = 100

    def dummy_1(self):
        """Dummy function to please pylint.

        """

    def dummy_2(self):
        """Dummy function to please pylint.

        """


class BasicPhysicsParameters(BasicParameters):
    """Class for basic physics testing parameters.

    """
    def __init__(self):
        super().__init__()
        self.temperature_reference = 3e2
        self.temperature_scale = 1e2
        self.boltzmann_constant = 8.314462618


class BasicSingleChainParameters(BasicPhysicsParameters):
    """Class for basic single-chain testing parameters.

    """
    def __init__(self):
        super().__init__()
        self.number_of_links_minimum = 5
        self.number_of_links_maximum = 25
        self.link_length_reference = 1.0
        self.link_length_scale = 1.0
        self.hinge_mass_reference = 1.0
        self.hinge_mass_scale = 1.0


class ExtraSingleChainParameters(BasicSingleChainParameters):
    """Class for extra single-chain testing parameters.

    """
    def __init__(self):
        super().__init__()
        self.link_stiffness_reference = 5e5
        self.link_stiffness_scale = 99e4
        self.nondimensional_link_stiffness_large = 1e4
        self.well_width_reference = 99e-2
        self.well_width_scale = 5e-1


class Parameters(ExtraSingleChainParameters):
    """Class for testing parameters.

    """
    def __init__(self):
        super().__init__()
        self.nondimensional_force_reference = 5e1
        self.nondimensional_force_scale = 1e2
