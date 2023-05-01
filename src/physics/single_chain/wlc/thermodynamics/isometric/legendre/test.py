"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters

parameters = Parameters()
WLC = physics.single_chain.wlc.thermodynamics.isometric.legendre.WLC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = WLC(
                parameters.number_of_links_minimum,
                parameters.link_length_reference,
                parameters.hinge_mass_reference,
                parameters.persistance_length_reference
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
                WLC(
                    number_of_links,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.persistance_length_reference
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
                WLC(
                    parameters.number_of_links_minimum,
                    link_length,
                    parameters.hinge_mass_reference,
                    parameters.persistance_length_reference
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
                WLC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass,
                    parameters.persistance_length_reference
                ).hinge_mass
            )

    def test_persistance_length(self):
        """Function to test the persistance length during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            self.assertEqual(
                persistance_length,
                WLC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    persistance_length
                ).persistance_length
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
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
                persistance_length,
                model.persistance_length
            )


class Nondimensional(unittest.TestCase):
    """Class for nondimensionalization tests.

    """
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            gibbs_free_energy = \
                model.gibbs_free_energy(
                    np.array(end_to_end_length),
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_gibbs_free_energy_per_link = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            gibbs_free_energy_per_link = \
                model.gibbs_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy_per_link / \
                parameters.boltzmann_constant/temperature \
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_relative_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            relative_gibbs_free_energy = \
                model.relative_gibbs_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_gibbs_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_relative_gibbs_free_energy
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_relative_gibbs_free_energy_per_link = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            relative_gibbs_free_energy_per_link = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_gibbs_free_energy_per_link / \
                parameters.boltzmann_constant/temperature \
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


class PerLink(unittest.TestCase):
    """Class for per-linkness tests.

    """
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            gibbs_free_energy = \
                model.gibbs_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            gibbs_free_energy_per_link = \
                model.gibbs_free_energy_per_link(
                    np.array(end_to_end_length),
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            relative_gibbs_free_energy = \
                model.relative_gibbs_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_gibbs_free_energy_per_link = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(end_to_end_length),
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            nondimensional_relative_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_gibbs_free_energy_per_link = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            gibbs_free_energy = \
                model.gibbs_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            gibbs_free_energy_0 = \
                model.gibbs_free_energy(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            relative_gibbs_free_energy = \
                model.relative_gibbs_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy \
                - gibbs_free_energy_0 \
                - relative_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                relative_gibbs_free_energy
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            gibbs_free_energy_per_link = \
                model.gibbs_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            gibbs_free_energy_per_link_0 = \
                model.gibbs_free_energy_per_link(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            relative_gibbs_free_energy_per_link = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                gibbs_free_energy_per_link \
                - gibbs_free_energy_per_link_0 \
                - relative_gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                relative_gibbs_free_energy_per_link
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_gibbs_free_energy = \
                model.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_gibbs_free_energy_0 = \
                model.nondimensional_gibbs_free_energy(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_gibbs_free_energy \
                - nondimensional_gibbs_free_energy_0 \
                - nondimensional_relative_gibbs_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_relative_gibbs_free_energy
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_end_to_end_length_per_link = \
                parameters. \
                nondimensional_end_to_end_length_per_link_reference + \
                parameters. \
                nondimensional_end_to_end_length_per_link_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_gibbs_free_energy_per_link = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link_0 = \
                model.nondimensional_gibbs_free_energy_per_link(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_gibbs_free_energy_per_link = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_gibbs_free_energy_per_link \
                - nondimensional_gibbs_free_energy_per_link_0 \
                - nondimensional_relative_gibbs_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_relative_gibbs_free_energy_per_link
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_gibbs_free_energy_0 = \
                model.relative_gibbs_free_energy(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            self.assertLessEqual(
                np.abs(relative_gibbs_free_energy_0),
                parameters.boltzmann_constant*temperature *
                number_of_links*parameters.zero
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_gibbs_free_energy_per_link_0 = \
                model.relative_gibbs_free_energy_per_link(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            self.assertLessEqual(
                np.abs(relative_gibbs_free_energy_per_link_0),
                parameters.boltzmann_constant*temperature *
                parameters.zero
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the zero
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_relative_gibbs_free_energy_0 = \
                model.nondimensional_relative_gibbs_free_energy(
                    np.array(parameters.zero)
                )
            self.assertLessEqual(
                np.abs(nondimensional_relative_gibbs_free_energy_0),
                number_of_links*parameters.zero
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
            persistance_length = \
                parameters.persistance_length_reference + \
                parameters.persistance_length_scale*(0.5 - np.random.rand())
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_relative_gibbs_free_energy_per_link_0 = \
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(parameters.zero)
                )
            self.assertLessEqual(
                np.abs(
                    nondimensional_relative_gibbs_free_energy_per_link_0
                ), parameters.zero
            )
