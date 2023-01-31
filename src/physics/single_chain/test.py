"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics


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
        self.number_of_loops = 8
        self.one = 1e0
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
        self.log_log_tol = 5e-2
        self.log_log_scale = 12e-1
        self.rel_tol_thermodynamic_limit = 1e-1
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
        self.nondimensional_link_stiffness_medium = 1e1
        self.well_width_reference = 99e-2
        self.well_width_scale = 5e-1
        self.nondimensional_well_width_small = 1e-2


class Parameters(ExtraSingleChainParameters):
    """Class for testing parameters.

    """
    def __init__(self):
        super().__init__()
        self.nondimensional_end_to_end_length_per_link_reference = 5e-1
        self.nondimensional_end_to_end_length_per_link_scale = 99e-2
        self.nondimensional_end_to_end_length_per_link_small = 25e-2
        self.nondimensional_force_reference = 5e1
        self.nondimensional_force_scale = 1e2
        self.nondimensional_force_small = 75e-2


parameters = Parameters()
Ideal = physics.single_chain.ideal.Ideal
FJC = physics.single_chain.fjc.FJC
EFJC = physics.single_chain.efjc.EFJC
SWFJC = physics.single_chain.swfjc.SWFJC


class FjcIdeal(unittest.TestCase):
    """Class for FJC becoming Ideal tests.

    """
    def test_end_to_end_length(self):
        """Function to test the FJC -> Ideal limit
        of the end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_fjc = \
                fjc.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            end_to_end_length_ideal = \
                ideal.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_fjc \
                - end_to_end_length_ideal
            residual_rel = \
                residual_abs / \
                end_to_end_length_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            end_to_end_length_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_per_link_fjc \
                - end_to_end_length_per_link_ideal
            residual_rel = \
                residual_abs / \
                end_to_end_length_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_end_to_end_length_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            nondimensional_end_to_end_length_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length_fjc \
                - nondimensional_end_to_end_length_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_end_to_end_length_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_end_to_end_length_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length_per_link_fjc \
                - nondimensional_end_to_end_length_per_link_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_force(self):
        """Function to test the FJC -> Ideal limit
        of the force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = \
                nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            force_fjc = \
                fjc.thermodynamics.isometric. \
                force(
                    np.array(end_to_end_length),
                    temperature
                )
            force_ideal = \
                ideal.thermodynamics.isometric. \
                force(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                force_fjc \
                - force_ideal
            residual_rel = \
                residual_abs / \
                force_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_nondimensional_force(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_force_fjc = \
                fjc.thermodynamics.isometric. \
                nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_force_ideal = \
                ideal.thermodynamics.isometric. \
                nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_force_fjc \
                - nondimensional_force_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_force_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_helmholtz_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = \
                nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            helmholtz_free_energy_fjc = \
                fjc.thermodynamics.isometric. \
                helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_ideal = \
                ideal.thermodynamics.isometric. \
                helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy_fjc \
                - helmholtz_free_energy_ideal
            residual_rel = \
                residual_abs / \
                helmholtz_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_helmholtz_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = \
                nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            helmholtz_free_energy_per_link_fjc = \
                fjc.thermodynamics.isometric. \
                helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_per_link_ideal = \
                ideal.thermodynamics.isometric. \
                helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy_per_link_fjc \
                - helmholtz_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                helmholtz_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = \
                nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            relative_helmholtz_free_energy_fjc = \
                fjc.thermodynamics.isometric. \
                relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy_ideal = \
                ideal.thermodynamics.isometric. \
                relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_helmholtz_free_energy_fjc \
                - relative_helmholtz_free_energy_ideal
            residual_rel = \
                residual_abs / \
                relative_helmholtz_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_relative_helmholtz_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = \
                nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            relative_helmholtz_free_energy_per_link_fjc = \
                fjc.thermodynamics.isometric. \
                relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy_per_link_ideal = \
                ideal.thermodynamics.isometric. \
                relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_helmholtz_free_energy_per_link_fjc \
                - relative_helmholtz_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                relative_helmholtz_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_nondimensional_helmholtz_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy_fjc = \
                fjc.thermodynamics.isometric. \
                nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_helmholtz_free_energy_ideal = \
                ideal.thermodynamics.isometric. \
                nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy_fjc \
                - nondimensional_helmholtz_free_energy_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_nondimensional_helmholtz_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy_per_link_fjc = \
                fjc.thermodynamics.isometric. \
                nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_helmholtz_free_energy_per_link_ideal = \
                ideal.thermodynamics.isometric. \
                nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy_per_link_fjc \
                - nondimensional_helmholtz_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_relative_helmholtz_free_energy_fjc = \
                fjc.thermodynamics.isometric. \
                nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy_ideal = \
                ideal.thermodynamics.isometric. \
                nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_relative_helmholtz_free_energy_fjc \
                - nondimensional_relative_helmholtz_free_energy_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_relative_helmholtz_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_nondimensional_relative_helmholtz_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_relative_helmholtz_free_energy_per_link_fjc = \
                fjc.thermodynamics.isometric. \
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy_per_link_ideal = \
                ideal.thermodynamics.isometric. \
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_relative_helmholtz_free_energy_per_link_fjc \
                - nondimensional_relative_helmholtz_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_relative_helmholtz_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_gibbs_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = \
                nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy_fjc = \
                fjc.thermodynamics.isotensional. \
                gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_ideal = \
                ideal.thermodynamics.isotensional. \
                gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy_fjc \
                - gibbs_free_energy_ideal
            residual_rel = \
                residual_abs / \
                gibbs_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_gibbs_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = \
                nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy_per_link_fjc \
                - gibbs_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                gibbs_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_relative_gibbs_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = \
                nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            relative_gibbs_free_energy_fjc = \
                fjc.thermodynamics.isotensional. \
                relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_ideal = \
                ideal.thermodynamics.isotensional. \
                relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                relative_gibbs_free_energy_fjc \
                - relative_gibbs_free_energy_ideal
            residual_rel = \
                residual_abs / \
                relative_gibbs_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_relative_gibbs_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = \
                nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            relative_gibbs_free_energy_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                relative_gibbs_free_energy_per_link_fjc \
                - relative_gibbs_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                relative_gibbs_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_gibbs_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_gibbs_free_energy_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_gibbs_free_energy_fjc \
                - nondimensional_gibbs_free_energy_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_gibbs_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_gibbs_free_energy_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_gibbs_free_energy_per_link_fjc \
                - nondimensional_gibbs_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_relative_gibbs_free_energy_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force)
                )
            nondimensional_relative_gibbs_free_energy_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_relative_gibbs_free_energy_fjc \
                - nondimensional_relative_gibbs_free_energy_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_relative_gibbs_free_energy_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_relative_gibbs_free_energy_per_link(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_relative_gibbs_free_energy_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_relative_gibbs_free_energy_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_relative_gibbs_free_energy_per_link_fjc \
                - nondimensional_relative_gibbs_free_energy_per_link_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_relative_gibbs_free_energy_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_equilibrium_distribution(self):
        """Function to test the FJC -> Ideal limit
        of the equilibrium distribution.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            end_to_end_length = \
                nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            equilibrium_distribution_fjc = \
                fjc.thermodynamics.isometric. \
                equilibrium_distribution(
                    np.array(end_to_end_length)
                )
            equilibrium_distribution_ideal = \
                ideal.thermodynamics.isometric. \
                equilibrium_distribution(
                    np.array(end_to_end_length)
                )
            residual_abs = \
                equilibrium_distribution_fjc \
                - equilibrium_distribution_ideal
            residual_rel = \
                residual_abs / \
                equilibrium_distribution_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_nondimensional_equilibrium_distribution(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional equilibrium distribution.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_equilibrium_distribution_fjc = \
                fjc.thermodynamics.isometric. \
                nondimensional_equilibrium_distribution(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_equilibrium_distribution_ideal = \
                ideal.thermodynamics.isometric. \
                nondimensional_equilibrium_distribution(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_equilibrium_distribution_fjc \
                - nondimensional_equilibrium_distribution_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_equilibrium_distribution_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_equilibrium_radial_distribution(self):
        """Function to test the FJC -> Ideal limit
        of the equilibrium radial distribution.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            end_to_end_length = \
                nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            equilibrium_radial_distribution_fjc = \
                fjc.thermodynamics.isometric. \
                equilibrium_radial_distribution(
                    np.array(end_to_end_length)
                )
            equilibrium_radial_distribution_ideal = \
                ideal.thermodynamics.isometric. \
                equilibrium_radial_distribution(
                    np.array(end_to_end_length)
                )
            residual_abs = \
                equilibrium_radial_distribution_fjc \
                - equilibrium_radial_distribution_ideal
            residual_rel = \
                residual_abs / \
                equilibrium_radial_distribution_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )

    def test_nondimensional_equilibrium_radial_distribution(self):
        """Function to test the FJC -> Ideal limit
        of the nondimensional equilibrium radial distribution.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_end_to_end_length_per_link = \
                parameters.nondimensional_end_to_end_length_per_link_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_equilibrium_radial_distribution_fjc = \
                fjc.thermodynamics.isometric. \
                nondimensional_equilibrium_radial_distribution(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_equilibrium_radial_distribution_ideal = \
                ideal.thermodynamics.isometric. \
                nondimensional_equilibrium_radial_distribution(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_equilibrium_radial_distribution_fjc \
                - nondimensional_equilibrium_radial_distribution_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_equilibrium_radial_distribution_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_end_to_end_length_per_link
            )


class EfjcIdeal(unittest.TestCase):
    """Class for EFJC becoming Ideal tests.

    """
    def test_end_to_end_length(self):
        """Function to test the EFJC -> Ideal limit
        of the end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = parameters.link_stiffness_scale
            efjc = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_efjc = \
                efjc.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            end_to_end_length_ideal = \
                ideal.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_efjc \
                - end_to_end_length_ideal
            residual_rel = \
                residual_abs / \
                end_to_end_length_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the EFJC -> Ideal limit
        of the end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = parameters.link_stiffness_scale
            efjc = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link_efjc = \
                efjc.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            end_to_end_length_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_per_link_efjc \
                - end_to_end_length_per_link_ideal
            residual_rel = \
                residual_abs / \
                end_to_end_length_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the EFJC -> Ideal limit
        of the nondimensional end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = parameters.link_stiffness_scale
            efjc = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_end_to_end_length_efjc = \
                efjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_end_to_end_length_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length_efjc \
                - nondimensional_end_to_end_length_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the EFJC -> Ideal limit
        of the nondimensional end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = parameters.link_stiffness_scale
            efjc = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_end_to_end_length_per_link_efjc = \
                efjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_end_to_end_length_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length_per_link_efjc \
                - nondimensional_end_to_end_length_per_link_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )


class SwfjcIdeal(unittest.TestCase):
    """Class for SWFJC becoming Ideal tests.

    """
    def test_end_to_end_length(self):
        """Function to test the SWFJC -> Ideal limit
        of the end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_swfjc = \
                swfjc.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            end_to_end_length_ideal = \
                ideal.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_swfjc \
                - end_to_end_length_ideal
            residual_rel = \
                residual_abs / \
                end_to_end_length_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the SWFJC -> Ideal limit
        of the end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link_swfjc = \
                swfjc.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            end_to_end_length_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_per_link_swfjc \
                - end_to_end_length_per_link_ideal
            residual_rel = \
                residual_abs / \
                end_to_end_length_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the SWFJC -> Ideal limit
        of the nondimensional end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_end_to_end_length_swfjc = \
                swfjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            nondimensional_end_to_end_length_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length_swfjc \
                - nondimensional_end_to_end_length_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the SWFJC -> Ideal limit
        of the nondimensional end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            ideal = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_end_to_end_length_per_link_swfjc = \
                swfjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_end_to_end_length_per_link_ideal = \
                ideal.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length_per_link_swfjc \
                - nondimensional_end_to_end_length_per_link_ideal
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_per_link_ideal
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                nondimensional_force
            )


class SwfjcFjc(unittest.TestCase):
    """Class for FJC becoming SWFJC tests.

    """
    def test_end_to_end_length(self):
        """Function to test the FJC -> SWFJC limit
        of the end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_fjc = \
                fjc.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            end_to_end_length_swfjc = \
                swfjc.thermodynamics.isotensional. \
                end_to_end_length(
                    np.array(force),
                    temperature
                )
            residual = \
                end_to_end_length_swfjc \
                - end_to_end_length_fjc
            self.assertLessEqual(
                np.abs(residual),
                number_of_links*well_width
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the FJC -> SWFJC limit
        of the end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            end_to_end_length_per_link_swfjc = \
                swfjc.thermodynamics.isotensional. \
                end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            residual = \
                end_to_end_length_per_link_swfjc \
                - end_to_end_length_per_link_fjc
            self.assertLessEqual(
                np.abs(residual),
                well_width
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the FJC -> SWFJC limit
        of the nondimensional end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_end_to_end_length_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            nondimensional_end_to_end_length_swfjc = \
                swfjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            residual = \
                nondimensional_end_to_end_length_swfjc \
                - nondimensional_end_to_end_length_fjc
            self.assertLessEqual(
                np.abs(residual),
                number_of_links*parameters.nondimensional_well_width_small
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the FJC -> SWFJC limit
        of the nondimensional end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            well_width = parameters.nondimensional_well_width_small*link_length
            fjc = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            swfjc = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            nondimensional_force = \
                parameters.nondimensional_force_small * \
                (1.0 - 0.5*np.random.rand())
            nondimensional_end_to_end_length_per_link_fjc = \
                fjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_end_to_end_length_per_link_swfjc = \
                swfjc.thermodynamics.isotensional. \
                nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            residual = \
                nondimensional_end_to_end_length_per_link_swfjc \
                - nondimensional_end_to_end_length_per_link_fjc
            self.assertLessEqual(
                np.abs(residual),
                parameters.nondimensional_well_width_small
            )
