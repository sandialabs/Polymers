"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters
from .....test import integrate

parameters = Parameters()
LENNARDJONESFJC = physics.single_chain.ufjc.lennard_jones.thermodynamics.\
    isotensional.LENNARDJONESFJC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = LENNARDJONESFJC(
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
                LENNARDJONESFJC(
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
                LENNARDJONESFJC(
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
                LENNARDJONESFJC(
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
                LENNARDJONESFJC(
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
            model = LENNARDJONESFJC(
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
    def test_end_to_end_length(self):
        """Function to test the nondimensionalization
        of the end-to-end length.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length = \
                model.nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.end_to_end_length(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length / \
                link_length \
                - nondimensional_end_to_end_length
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the nondimensionalization
        of the end-to-end length per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length_per_link = \
                model.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link = \
                model.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_per_link / \
                link_length \
                - nondimensional_end_to_end_length_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_gibbs_free_energy(self):
        """Function to test the nondimensionalization
        of the Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy = \
                model.gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_gibbs_free_energy_per_link(self):
        """Function to test the nondimensionalization
        of the Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy = \
                model.gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_relative_gibbs_free_energy(self):
        """Function to test the nondimensionalization
        of the relative Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy = \
                model.relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_relative_gibbs_free_energy_per_link(self):
        """Function to test the nondimensionalization
        of the relative Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )


class PerLink(unittest.TestCase):
    """Class for per-linkness tests.

    """
    def test_end_to_end_length(self):
        """Function to test the per-linkness
        of the end-to-end length.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.end_to_end_length(
                    np.array(force),
                    temperature
                )
            end_to_end_length_per_link = \
                model.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length / \
                number_of_links \
                - end_to_end_length_per_link
            residual_rel = \
                residual_abs / \
                end_to_end_length_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the per-linkness
        of the nondimensional end-to-end length.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length = \
                model.nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_end_to_end_length_per_link = \
                model.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_end_to_end_length / \
                number_of_links \
                - nondimensional_end_to_end_length_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_end_to_end_length_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_gibbs_free_energy(self):
        """Function to test the per-linkness
        of the Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy = \
                model.gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_per_link = \
                model.gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy / \
                number_of_links \
                - gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_relative_gibbs_free_energy(self):
        """Function to test the per-linkness
        of the relative Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            relative_gibbs_free_energy = \
                model.relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_per_link = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                relative_gibbs_free_energy / \
                number_of_links \
                - relative_gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                relative_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_gibbs_free_energy(self):
        """Function to test the per-linkness
        of the nondimensional Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_gibbs_free_energy / \
                number_of_links \
                - nondimensional_gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the per-linkness
        of the nondimensional relative Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_relative_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy_per_link = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_relative_gibbs_free_energy / \
                number_of_links \
                - nondimensional_relative_gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_relative_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )


class Relative(unittest.TestCase):
    """Class for relativeness tests.

    """
    def test_gibbs_free_energy(self):
        """Function to test the relativeness
        of the Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy = \
                model.gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_0 = \
                model.gibbs_free_energy(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            relative_gibbs_free_energy = \
                model.relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy \
                - gibbs_free_energy_0 \
                - relative_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                gibbs_free_energy_0
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature/link_length
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_gibbs_free_energy_per_link(self):
        """Function to test the relativeness
        of the Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            gibbs_free_energy_per_link = \
                model.gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_per_link_0 = \
                model.gibbs_free_energy_per_link(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            relative_gibbs_free_energy_per_link = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy_per_link \
                - gibbs_free_energy_per_link_0 \
                - relative_gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                gibbs_free_energy_per_link_0
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature/link_length
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_gibbs_free_energy(self):
        """Function to test the relativeness
        of the nondimensional Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_0 = \
                model.nondimensional_gibbs_free_energy(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_gibbs_free_energy \
                - nondimensional_gibbs_free_energy_0 \
                - nondimensional_relative_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy_0
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_gibbs_free_energy_per_link(self):
        """Function to test the relativeness
        of the nondimensional Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_gibbs_free_energy_per_link = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link_0 = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy_per_link = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            residual_abs = \
                nondimensional_gibbs_free_energy_per_link \
                - nondimensional_gibbs_free_energy_per_link_0 \
                - nondimensional_relative_gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy_per_link_0
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )


class Zero(unittest.TestCase):
    """Class for zero tests.

    """
    def test_relative_gibbs_free_energy(self):
        """Function to test the zero
        of the relative Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_gibbs_free_energy_0 = \
                model.relative_gibbs_free_energy(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            self.assertLessEqual(
                relative_gibbs_free_energy_0,
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature
            )

    def test_relative_gibbs_free_energy_per_link(self):
        """Function to test the zero
        of the relative Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_gibbs_free_energy_per_link_0 = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(
                        parameters.zero *
                        parameters.boltzmann_constant*temperature/link_length
                    ),
                    temperature
                )
            self.assertLessEqual(
                relative_gibbs_free_energy_per_link_0,
                parameters.abs_tol *
                parameters.boltzmann_constant*temperature
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the zero
        of the relative Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_relative_gibbs_free_energy_0 = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(parameters.zero),
                    temperature
                )
            self.assertLessEqual(
                nondimensional_relative_gibbs_free_energy_0,
                parameters.abs_tol
            )

    def test_nondimensional_relative_gibbs_free_energy_per_link(self):
        """Function to test the zero
        of the nondimensional relative Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_relative_gibbs_free_energy_per_link_0 = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(parameters.zero),
                    temperature
                )
            self.assertLessEqual(
                nondimensional_relative_gibbs_free_energy_per_link_0,
                parameters.abs_tol
            )


class Connection(unittest.TestCase):
    """Class for connection tests.

    """
    def test_end_to_end_length(self):
        """Function to test the connection
        of the end-to-end length.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.end_to_end_length(
                    np.array(force),
                    temperature
                )
            h_step = parameters.rel_tol * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_from_derivative = -(
                model.relative_gibbs_free_energy(
                    np.array(force + 0.5*h_step),
                    temperature
                )
                - model.relative_gibbs_free_energy(
                    np.array(force - 0.5*h_step),
                    temperature
                ))/h_step
            residual_abs = \
                end_to_end_length \
                - end_to_end_length_from_derivative
            self.assertLessEqual(
                np.abs(residual_abs), h_step
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the connection
        of the end-to-end length per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link = \
                model.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            h_step = parameters.rel_tol * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link_from_derivative = -(
                model.relative_gibbs_free_energy_per_link(
                    np.array(force + 0.5*h_step),
                    temperature
                )
                - model.relative_gibbs_free_energy_per_link(
                    np.array(force - 0.5*h_step),
                    temperature
                ))/h_step
            residual_abs = \
                end_to_end_length_per_link \
                - end_to_end_length_per_link_from_derivative
            self.assertLessEqual(
                np.abs(residual_abs), h_step
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the connection
        of the nondimensional end-to-end length.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length = \
                model.nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                    temperature
                )
            h_step = parameters.rel_tol
            nondimensional_end_to_end_length_from_derivative = -(
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force + 0.5*h_step),
                    temperature
                )
                - model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force - 0.5*h_step),
                    temperature
                ))/h_step
            residual_abs = \
                nondimensional_end_to_end_length \
                - nondimensional_end_to_end_length_from_derivative
            self.assertLessEqual(
                np.abs(residual_abs), h_step
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the connection
        of the nondimensional end-to-end length per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length_per_link = \
                model.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            h_step = parameters.rel_tol
            nondimensional_end_to_end_length_per_link_from_derivative = -(
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force + 0.5*h_step),
                    temperature
                )
                - model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force - 0.5*h_step),
                    temperature
                ))/h_step
            residual_abs = \
                nondimensional_end_to_end_length_per_link \
                - nondimensional_end_to_end_length_per_link_from_derivative
            self.assertLessEqual(
                np.abs(residual_abs), h_step
            )


class Legendre(unittest.TestCase):
    """Class for Legendre tests.

    """
    def test_gibbs_free_energy(self):
        """Function to test the Legendre transformation
        of the Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.end_to_end_length(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy = \
                model.gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_legendre = \
                model.legendre.helmholtz_free_energy(
                    np.array(force),
                    temperature
                ) - force*end_to_end_length
            residual_abs = \
                gibbs_free_energy \
                - gibbs_free_energy_legendre
            residual_rel = \
                residual_abs / \
                gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_gibbs_free_energy_per_link(self):
        """Function to test the Legendre transformation
        of the Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link = \
                model.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_per_link = \
                model.gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_per_link_legendre = \
                model.legendre.helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                ) - force*end_to_end_length_per_link
            residual_abs = \
                gibbs_free_energy_per_link \
                - gibbs_free_energy_per_link_legendre
            residual_rel = \
                residual_abs / \
                gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_relative_gibbs_free_energy(self):
        """Function to test the Legendre transformation
        of the relative Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.end_to_end_length(
                    np.array(force),
                    temperature
                )
            end_to_end_length_0 = \
                model.end_to_end_length(
                    np.array(
                        parameters.boltzmann_constant*temperature/link_length
                        * parameters.zero
                    ),
                    temperature
                )
            relative_gibbs_free_energy = \
                model.relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_legendre = \
                model.legendre.relative_helmholtz_free_energy(
                    np.array(force),
                    temperature
                ) - force*end_to_end_length \
                + end_to_end_length_0 \
                * parameters.boltzmann_constant*temperature/link_length \
                * parameters.zero
            residual_abs = \
                relative_gibbs_free_energy \
                - relative_gibbs_free_energy_legendre
            residual_rel = \
                residual_abs / \
                relative_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_relative_gibbs_free_energy_per_link(self):
        """Function to test the Legendre transformation
        of the relative Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length_per_link = \
                model.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            end_to_end_length_per_link_0 = \
                model.end_to_end_length_per_link(
                    np.array(
                        parameters.boltzmann_constant*temperature/link_length
                        * parameters.zero
                    ),
                    temperature
                )
            relative_gibbs_free_energy_per_link = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_per_link_legendre = \
                model.legendre.relative_helmholtz_free_energy_per_link(
                    np.array(force),
                    temperature
                ) - force*end_to_end_length_per_link \
                + end_to_end_length_per_link_0 \
                * parameters.boltzmann_constant*temperature/link_length \
                * parameters.zero
            residual_abs = \
                relative_gibbs_free_energy_per_link \
                - relative_gibbs_free_energy_per_link_legendre
            residual_rel = \
                residual_abs / \
                relative_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_gibbs_free_energy(self):
        """Function to test the Legendre transformation
        of the nondimensional Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length = \
                model.nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_legendre = \
                model.legendre.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                ) - nondimensional_force*nondimensional_end_to_end_length
            residual_abs = \
                nondimensional_gibbs_free_energy \
                - nondimensional_gibbs_free_energy_legendre
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_gibbs_free_energy_per_link(self):
        """Function to test the Legendre transformation
        of the nondimensional Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length_per_link = \
                model.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link_legendre = \
                model.legendre.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                ) - nondimensional_force * \
                nondimensional_end_to_end_length_per_link
            residual_abs = \
                nondimensional_gibbs_free_energy_per_link \
                - nondimensional_gibbs_free_energy_per_link_legendre
            residual_rel = \
                residual_abs / \
                nondimensional_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the Legendre transformation
        of the nondimensional relative Gibbs free energy.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length = \
                model.nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_end_to_end_length_0 = \
                model.nondimensional_end_to_end_length(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy_legendre = \
                model.legendre.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_force),
                    temperature
                ) - nondimensional_force*nondimensional_end_to_end_length \
                + parameters.zero*nondimensional_end_to_end_length_0
            residual_abs = \
                nondimensional_relative_gibbs_free_energy \
                - nondimensional_relative_gibbs_free_energy_legendre
            residual_rel = \
                residual_abs / \
                nondimensional_relative_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_relative_gibbs_free_energy_per_link(self):
        """Function to test the Legendre transformation
        of the nondimensional relative Gibbs free energy per link.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            nondimensional_end_to_end_length_per_link = \
                model.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_end_to_end_length_per_link_0 = \
                model.nondimensional_end_to_end_length_per_link(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy_per_link = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy_per_link_legendre = \
                model.legendre.\
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                ) - nondimensional_force * \
                nondimensional_end_to_end_length_per_link \
                + parameters.zero*nondimensional_end_to_end_length_per_link_0
            residual_abs = \
                nondimensional_relative_gibbs_free_energy_per_link \
                - nondimensional_relative_gibbs_free_energy_per_link_legendre
            residual_rel = \
                residual_abs / \
                nondimensional_relative_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )


class LegendreConnection(unittest.TestCase):
    """Class for Legendre connection tests.

    """
    def test_force(self):
        """Function to test the Legendre transformation connection
        of the force.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            h_step = parameters.rel_tol * \
                parameters.boltzmann_constant*temperature/link_length
            force_from_derivative = (
                model.legendre.relative_helmholtz_free_energy(
                    np.array(
                        force + 0.5*h_step
                    ), temperature
                )
                - model.legendre.relative_helmholtz_free_energy(
                    np.array(
                        force - 0.5*h_step
                    ), temperature
                ))/(
                model.end_to_end_length(
                    np.array(
                        force + 0.5*h_step
                    ), temperature
                )
                - model.end_to_end_length(
                    np.array(
                        force - 0.5*h_step
                    ), temperature
                ))
            residual_abs = \
                force \
                - force_from_derivative
            residual_rel = residual_abs/force
            self.assertLessEqual(
                np.abs(residual_rel), h_step
            )

    def test_nondimensional_force(self):
        """Function to test the Legendre transformation connection
        of the nondimensional force.

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
            model = LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            lambda_max = (13/7)**(1/6)
            nondimensional_force_max = link_stiffness / \
                parameters.boltzmann_constant/temperature * \
                link_length**2/6.0*(1/lambda_max**7 - 1/lambda_max**13)
            nondimensional_force = nondimensional_force_max*np.random.rand()
            h_step = parameters.rel_tol
            nondimensional_force_from_derivative = (
                model.legendre.
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(
                        nondimensional_force + 0.5*h_step
                    ), temperature
                )
                - model.legendre.
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(
                        nondimensional_force - 0.5*h_step
                    ), temperature
                ))/(
                model.nondimensional_end_to_end_length_per_link(
                    np.array(
                        nondimensional_force + 0.5*h_step
                    ), temperature
                )
                - model.nondimensional_end_to_end_length_per_link(
                    np.array(
                        nondimensional_force - 0.5*h_step
                    ), temperature
                ))
            residual_abs = \
                nondimensional_force \
                - nondimensional_force_from_derivative
            residual_rel = residual_abs/nondimensional_force
            self.assertLessEqual(
                np.abs(residual_rel), h_step
            )


class Asymptotic(unittest.TestCase):
    """Class for asymptotic tests.

    """
    def test_end_to_end_length(self):
        """Function to test the asymptotics
        of the end-to-end length.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length(
                            force, temperature
                        ) -
                        model.asymptotic.end_to_end_length(
                            force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length(
                            force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(0.5*log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the asymptotics
        of the end-to-end length per link.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length_per_link(
                            force, temperature
                        ) -
                        model.asymptotic.end_to_end_length_per_link(
                            force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length_per_link(
                            force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(0.5*log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the asymptotics
        of the nondimensional end-to-end length.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length(
                            nondimensional_force, temperature
                        ) -
                        model.asymptotic.
                        nondimensional_end_to_end_length(
                            nondimensional_force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length(
                            nondimensional_force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(0.5*log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the asymptotics
        of the nondimensional end-to-end length per link.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length_per_link(
                            nondimensional_force, temperature
                        ) -
                        model.asymptotic.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length_per_link(
                            nondimensional_force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(0.5*log_log_slope + 1.0),
                parameters.log_log_tol
            )


class AsymptoticReduced(unittest.TestCase):
    """Class for asymptotic tests for reduced approach.

    """
    def test_end_to_end_length(self):
        """Function to test the asymptotics
        of the end-to-end length.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length(
                            force, temperature
                        ) -
                        model.asymptotic.reduced.end_to_end_length(
                            force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length(
                            force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                5*parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the asymptotics
        of the end-to-end length per link.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length_per_link(
                            force, temperature
                        ) -
                        model.asymptotic.reduced.end_to_end_length_per_link(
                            force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    force = nondimensional_force * \
                        parameters.boltzmann_constant*temperature/link_length
                    return (
                        model.end_to_end_length_per_link(
                            force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the asymptotics
        of the nondimensional end-to-end length.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length(
                            nondimensional_force, temperature
                        ) -
                        model.asymptotic.reduced.
                        nondimensional_end_to_end_length(
                            nondimensional_force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length(
                            nondimensional_force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the asymptotics
        of the nondimensional end-to-end length per link.

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
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_link_stiffness):
                link_stiffness = nondimensional_link_stiffness * \
                    parameters.boltzmann_constant*temperature/link_length**2
                model = LENNARDJONESFJC(
                    number_of_links,
                    link_length,
                    hinge_mass,
                    link_stiffness
                )

                def integrand_numerator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length_per_link(
                            nondimensional_force, temperature
                        ) -
                        model.asymptotic.reduced.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_force, temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_force):
                    return (
                        model.nondimensional_end_to_end_length_per_link(
                            nondimensional_force, temperature
                        )
                    )**2

                nondimensional_link_stretch_max = (13/7)**(1/6)
                nondimensional_force_max = nondimensional_link_stretch_max
                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    nondimensional_force_max,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_link_stiffness_big
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_link_stiffness_big *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )
