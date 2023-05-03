"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters

parameters = Parameters()
WLC = physics.single_chain.wlc.thermodynamics.WLC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = WLC(
                parameters.number_of_links_maximum,
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
                    parameters.number_of_links_maximum,
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
                    parameters.number_of_links_maximum,
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
                    parameters.number_of_links_maximum,
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
                model.hinge_mass,
                persistance_length
            )
            self.assertEqual(
                persistance_length,
                model.persistance_length
            )


class ThermodynamicLimit(unittest.TestCase):
    """Class for thermodynamic limit tests.

    """
    def test_end_to_end_length(self):
        """Function to test the thermodynamic limit
        of the end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            force = \
                model.isometric.force(
                    np.array(end_to_end_length),
                    temperature
                )
            end_to_end_length_out = \
                model.isotensional.end_to_end_length(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length \
                - end_to_end_length_out
            residual_rel = residual_abs/end_to_end_length
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the thermodynamic limit
        of the end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            force = \
                model.isometric.force(
                    np.array(end_to_end_length),
                    temperature
                )
            end_to_end_length_per_link = \
                nondimensional_end_to_end_length_per_link*link_length
            end_to_end_length_per_link_out = \
                model.isotensional.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            residual_abs = \
                end_to_end_length_per_link \
                - end_to_end_length_per_link_out
            residual_rel = residual_abs/end_to_end_length_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the thermodynamic limit
        of the nondimensional end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            nondimensional_force = \
                model.isometric.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_end_to_end_length = \
                nondimensional_end_to_end_length_per_link*number_of_links
            nondimensional_end_to_end_length_out = \
                model.isotensional.nondimensional_end_to_end_length(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length \
                - nondimensional_end_to_end_length_out
            residual_rel = residual_abs/nondimensional_end_to_end_length
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the thermodynamic limit
        of the nondimensional end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            nondimensional_force = \
                model.isometric.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_end_to_end_length_per_link_out = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            residual_abs = \
                nondimensional_end_to_end_length_per_link \
                - nondimensional_end_to_end_length_per_link_out
            residual_rel = residual_abs / \
                nondimensional_end_to_end_length_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_force(self):
        """Function to test the thermodynamic limit
        of the force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.isotensional.end_to_end_length(
                    np.array(force),
                    temperature
                )
            force_out = \
                model.isometric.force(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                force \
                - force_out
            residual_rel = residual_abs/force
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_force(self):
        """Function to test the thermodynamic limit
        of the nondimensional force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_force_out = \
                model.isometric.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_force \
                - nondimensional_force_out
            residual_rel = residual_abs/nondimensional_force
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_helmholtz_free_energy(self):
        """Function to test the thermodynamic limit
        of the Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            force = \
                model.isometric.force(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy = \
                model.isometric.helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_out = \
                model.isotensional.gibbs_free_energy(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length
            residual_abs = \
                helmholtz_free_energy \
                - helmholtz_free_energy_out \
                - parameters.boltzmann_constant*temperature*np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )
            residual_rel = residual_abs/helmholtz_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_helmholtz_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            force = \
                model.isometric.force(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_per_link = \
                model.isometric.helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_per_link_out = \
                model.isotensional.gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length/number_of_links
            residual_abs = \
                helmholtz_free_energy_per_link \
                - helmholtz_free_energy_per_link_out \
                - parameters.boltzmann_constant*temperature*np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )/number_of_links
            residual_rel = residual_abs/helmholtz_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the thermodynamic limit
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            force = \
                model.isometric.force(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy = \
                model.isometric.relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy_out = \
                model.isotensional.relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length
            residual_abs = \
                relative_helmholtz_free_energy \
                - relative_helmholtz_free_energy_out
            residual_rel = residual_abs/relative_helmholtz_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_relative_helmholtz_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            force = \
                model.isometric.force(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy_per_link = \
                model.isometric.relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy_per_link_out = \
                model.isotensional.relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length/number_of_links
            residual_abs = \
                relative_helmholtz_free_energy_per_link \
                - relative_helmholtz_free_energy_per_link_out
            residual_rel = residual_abs/relative_helmholtz_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_helmholtz_free_energy(self):
        """Function to test the thermodynamic limit
        of the nondimensional Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            nondimensional_force = \
                model.isometric.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_helmholtz_free_energy = \
                model.isometric.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_helmholtz_free_energy_out = \
                model.isotensional.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                ) + nondimensional_force * \
                nondimensional_end_to_end_length_per_link*number_of_links
            residual_abs = \
                nondimensional_helmholtz_free_energy \
                - nondimensional_helmholtz_free_energy_out \
                - np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )
            residual_rel = residual_abs / \
                nondimensional_helmholtz_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_helmholtz_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the nondimensional Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            nondimensional_force = \
                model.isometric.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_helmholtz_free_energy_per_link = \
                model.isometric.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_helmholtz_free_energy_per_link_out = \
                model.isotensional.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                ) + nondimensional_force * \
                nondimensional_end_to_end_length_per_link
            residual_abs = \
                nondimensional_helmholtz_free_energy_per_link \
                - nondimensional_helmholtz_free_energy_per_link_out \
                - np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )/number_of_links
            residual_rel = residual_abs / \
                nondimensional_helmholtz_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the thermodynamic limit
        of the nondimensional relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            nondimensional_force = \
                model.isometric.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy = \
                model.isometric. \
                nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy_out = \
                model.isotensional. \
                nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force)
                ) + nondimensional_force * \
                nondimensional_end_to_end_length_per_link*number_of_links
            residual_abs = \
                nondimensional_relative_helmholtz_free_energy \
                - nondimensional_relative_helmholtz_free_energy_out
            residual_rel = residual_abs / \
                nondimensional_relative_helmholtz_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_relative_helmholtz_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the nondimensional relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
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
            nondimensional_force = \
                model.isometric.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy_per_link = \
                model.isometric. \
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy_per_link_out = \
                model.isotensional. \
                nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force)
                ) + nondimensional_force * \
                nondimensional_end_to_end_length_per_link
            residual_abs = \
                nondimensional_relative_helmholtz_free_energy_per_link \
                - nondimensional_relative_helmholtz_free_energy_per_link_out
            residual_rel = residual_abs / \
                nondimensional_relative_helmholtz_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_gibbs_free_energy(self):
        """Function to test the thermodynamic limit
        of the Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.isotensional.end_to_end_length(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy = \
                model.isotensional.gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_out = \
                model.isometric.helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                ) - force*end_to_end_length
            residual_abs = \
                gibbs_free_energy \
                - gibbs_free_energy_out \
                + parameters.boltzmann_constant*temperature*np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )
            residual_rel = residual_abs/gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_gibbs_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.isotensional.end_to_end_length(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_per_link = \
                model.isotensional.gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            gibbs_free_energy_per_link_out = \
                model.isometric.helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                ) - force*end_to_end_length/number_of_links
            residual_abs = \
                gibbs_free_energy_per_link \
                - gibbs_free_energy_per_link_out \
                + parameters.boltzmann_constant*temperature*np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )/number_of_links
            residual_rel = residual_abs/gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_relative_gibbs_free_energy(self):
        """Function to test the thermodynamic limit
        of the relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.isotensional.end_to_end_length(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy = \
                model.isotensional.relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_out = \
                model.isometric.relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                ) - force*end_to_end_length
            residual_abs = \
                relative_gibbs_free_energy \
                - relative_gibbs_free_energy_out
            residual_rel = residual_abs/relative_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_relative_gibbs_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            force = nondimensional_force * \
                parameters.boltzmann_constant*temperature/link_length
            end_to_end_length = \
                model.isotensional.end_to_end_length(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_per_link = \
                model.isotensional.relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                )
            relative_gibbs_free_energy_per_link_out = \
                model.isometric.relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                ) - force*end_to_end_length/number_of_links
            residual_abs = \
                relative_gibbs_free_energy_per_link \
                - relative_gibbs_free_energy_per_link_out
            residual_rel = residual_abs/relative_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_gibbs_free_energy(self):
        """Function to test the thermodynamic limit
        of the relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_gibbs_free_energy = \
                model.isotensional.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_out = \
                model.isometric.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                ) - nondimensional_force * \
                nondimensional_end_to_end_length_per_link*number_of_links
            residual_abs = \
                nondimensional_gibbs_free_energy \
                - nondimensional_gibbs_free_energy_out \
                + np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )
            residual_rel = residual_abs / \
                nondimensional_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_gibbs_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_gibbs_free_energy_per_link = \
                model.isotensional.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                )
            nondimensional_gibbs_free_energy_per_link_out = \
                model.isometric.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                ) - nondimensional_force * \
                nondimensional_end_to_end_length_per_link
            residual_abs = \
                nondimensional_gibbs_free_energy_per_link \
                - nondimensional_gibbs_free_energy_per_link_out \
                + np.log(
                    4.0*np.sin(np.arccos(np.exp(
                        -number_of_links*link_length/persistance_length
                    )))*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )/number_of_links
            residual_rel = residual_abs / \
                nondimensional_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the thermodynamic limit
        of the relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_relative_gibbs_free_energy = \
                model.isotensional. \
                nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force)
                )
            nondimensional_relative_gibbs_free_energy_out = \
                model.isometric. \
                nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                ) - nondimensional_force * \
                nondimensional_end_to_end_length_per_link*number_of_links
            residual_abs = \
                nondimensional_relative_gibbs_free_energy \
                - nondimensional_relative_gibbs_free_energy_out
            residual_rel = residual_abs / \
                nondimensional_relative_gibbs_free_energy
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )

    def test_nondimensional_relative_gibbs_free_energy_per_link(self):
        """Function to test the thermodynamic limit
        of the relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            persistance_length = \
                parameters.nondimensional_persistance_length_small * \
                number_of_links*link_length
            model = WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_relative_gibbs_free_energy_per_link = \
                model.isotensional. \
                nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force)
                )
            nondimensional_relative_gibbs_free_energy_per_link_out = \
                model.isometric. \
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                ) - nondimensional_force * \
                nondimensional_end_to_end_length_per_link
            residual_abs = \
                nondimensional_relative_gibbs_free_energy_per_link \
                - nondimensional_relative_gibbs_free_energy_per_link_out
            residual_rel = residual_abs / \
                nondimensional_relative_gibbs_free_energy_per_link
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol_thermodynamic_limit
            )
