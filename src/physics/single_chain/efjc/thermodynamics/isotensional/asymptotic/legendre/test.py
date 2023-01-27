"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters

parameters = Parameters()
EFJC = physics.single_chain.efjc.thermodynamics. \
    isotensional.asymptotic.legendre.EFJC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = EFJC(
                parameters.number_of_links_minimum,
                parameters.link_length_reference,
                parameters.hinge_mass_reference,
                parameters.link_stiffness_reference
            )

    def test_number_of_links(self):
        """Function to test the number of links during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            self.assertEqual(
                number_of_links,
                EFJC(
                    number_of_links,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.link_stiffness_reference
                ).number_of_links
            )

    def test_link_length(self):
        """Function to test the link length during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            self.assertEqual(
                link_length,
                EFJC(
                    parameters.number_of_links_minimum,
                    link_length,
                    parameters.hinge_mass_reference,
                    parameters.link_stiffness_reference
                ).link_length
            )

    def test_hinge_mass(self):
        """Function to test the hinge mass during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            self.assertEqual(
                hinge_mass,
                EFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass,
                    parameters.link_stiffness_reference
                ).hinge_mass
            )

    def test_link_stiffness(self):
        """Function to test the link stiffness during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            self.assertEqual(
                link_stiffness,
                EFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    link_stiffness
                ).link_stiffness
            )

    def test_all_parameters(self):
        """Function to test all parameters during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            self.assertEqual(
                number_of_links,
                model.number_of_links
            )
            self.assertEqual(
                link_length,
                model.link_length
            )
            self.assertEqual(
                hinge_mass,
                model.hinge_mass
            )
            self.assertEqual(
                link_stiffness,
                model.link_stiffness
            )


class Nondimensional(unittest.TestCase):
    """Class for nondimensionalization tests.

    """
    def test_helmholtz_free_energy(self):
        """Function to test the nondimensionalization
        of the Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            helmholtz_free_energy = \
                model.helmholtz_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_helmholtz_free_energy_per_link(self):
        """Function to test the nondimensionalization
        of the Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            helmholtz_free_energy = \
                model.helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the nondimensionalization
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            helmholtz_free_energy = \
                model.relative_helmholtz_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_relative_helmholtz_free_energy_per_link(self):
        """Function to test the nondimensionalization
        of the relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            helmholtz_free_energy = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )


class PerLink(unittest.TestCase):
    """Class for per-linkness tests.

    """
    def test_helmholtz_free_energy(self):
        """Function to test the per-linkness
        of the Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            helmholtz_free_energy = \
                model.helmholtz_free_energy(
                    np.array(force),
                    temperature
                )
            helmholtz_free_energy_per_link = \
                model.helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy / \
                number_of_links \
                - helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                helmholtz_free_energy_per_link
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the per-linkness
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            relative_helmholtz_free_energy = \
                model.relative_helmholtz_free_energy(
                    np.array(force),
                    temperature
                )
            relative_helmholtz_free_energy_per_link = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                relative_helmholtz_free_energy / \
                number_of_links \
                - relative_helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                relative_helmholtz_free_energy_per_link
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_nondimensional_helmholtz_free_energy(self):
        """Function to test the per-linkness
        of the nondimensional Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_helmholtz_free_energy_per_link = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy / \
                number_of_links \
                - nondimensional_helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy_per_link
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the per-linkness
        of the nondimensional relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_relative_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_relative_helmholtz_free_energy_per_link = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_relative_helmholtz_free_energy / \
                number_of_links \
                - nondimensional_relative_helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_relative_helmholtz_free_energy_per_link
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )


class Relative(unittest.TestCase):
    """Class for relativeness tests.

    """
    def test_helmholtz_free_energy(self):
        """Function to test the relativeness
        of the Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            helmholtz_free_energy = \
                model.helmholtz_free_energy(
                    np.array(force),
                    temperature
                )
            helmholtz_free_energy_0 = \
                model.helmholtz_free_energy(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            relative_helmholtz_free_energy = \
                model.relative_helmholtz_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy \
                - helmholtz_free_energy_0 \
                - relative_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                helmholtz_free_energy_0
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature/link_length
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_helmholtz_free_energy_per_link(self):
        """Function to test the relativeness
        of the Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            helmholtz_free_energy_per_link = \
                model.helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            helmholtz_free_energy_per_link_0 = \
                model.helmholtz_free_energy_per_link(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            relative_helmholtz_free_energy_per_link = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy_per_link \
                - helmholtz_free_energy_per_link_0 \
                - relative_helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                helmholtz_free_energy_per_link_0
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature/link_length
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_nondimensional_helmholtz_free_energy(self):
        """Function to test the relativeness
        of the nondimensional Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_helmholtz_free_energy_0 = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy \
                - nondimensional_helmholtz_free_energy_0 \
                - nondimensional_relative_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy_0
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

    def test_nondimensional_helmholtz_free_energy_per_link(self):
        """Function to test the relativeness
        of the nondimensional Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            nondimensional_force = \
                parameters.nondimensional_force_reference + \
                parameters.nondimensional_force_scale*(0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_helmholtz_free_energy_per_link = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_helmholtz_free_energy_per_link_0 = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_helmholtz_free_energy_per_link = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy_per_link \
                - nondimensional_helmholtz_free_energy_per_link_0 \
                - nondimensional_relative_helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_helmholtz_free_energy_per_link_0
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )


class Zero(unittest.TestCase):
    """Class for zero tests.

    """
    def test_relative_helmholtz_free_energy(self):
        """Function to test the zero
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_helmholtz_free_energy_0 = \
                model.relative_helmholtz_free_energy(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            self.assertLessEqual(
                relative_helmholtz_free_energy_0,
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature
            )

    def test_relative_helmholtz_free_energy_per_link(self):
        """Function to test the zero
        of the relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_helmholtz_free_energy_per_link_0 = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            self.assertLessEqual(
                relative_helmholtz_free_energy_per_link_0,
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the zero
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_relative_helmholtz_free_energy_0 = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(parameters.zero),
                    temperature
                )
            self.assertLessEqual(
                nondimensional_relative_helmholtz_free_energy_0,
                parameters.abs_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy_per_link(self):
        """Function to test the zero
        of the nondimensional relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                np.random.randint(
                    parameters.number_of_links_minimum,
                    high=parameters.number_of_links_maximum
                )
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            link_stiffness = \
                parameters.link_stiffness_reference + \
                parameters.link_stiffness_scale*(0.5 - np.random.rand())
            model = EFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_relative_helmholtz_free_energy_per_link_0 = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(parameters.zero),
                    temperature
                )
            self.assertLessEqual(
                nondimensional_relative_helmholtz_free_energy_per_link_0,
                parameters.abs_tol
            )
