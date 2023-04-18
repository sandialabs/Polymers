"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters

parameters = Parameters()
SWFJC = physics.single_chain.swfjc.thermodynamics.SWFJC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = SWFJC(
                parameters.number_of_links_minimum,
                parameters.link_length_reference,
                parameters.hinge_mass_reference,
                parameters.well_width_reference
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
                SWFJC(
                    number_of_links,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.well_width_reference
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
                SWFJC(
                    parameters.number_of_links_minimum,
                    link_length,
                    parameters.hinge_mass_reference,
                    parameters.well_width_reference
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
                SWFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass,
                    parameters.well_width_reference
                ).hinge_mass
            )

    def test_well_width(self):
        """Function to test the well width during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            self.assertEqual(
                well_width,
                SWFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    well_width
                ).well_width
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
                well_width,
                model.well_width
            )


class Legendre(unittest.TestCase):
    """Class for Legendre transformation tests.

    """
    def test_force(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
                model.isometric.legendre.force(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                force \
                - force_out
            residual_rel = residual_abs/force
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_force(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
                model.isometric.legendre.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_force \
                - nondimensional_force_out
            residual_rel = residual_abs/nondimensional_force
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_helmholtz_free_energy(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
            helmholtz_free_energy_legendre = \
                model.isotensional.gibbs_free_energy(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length
            helmholtz_free_energy_legendre_out = \
                model.isometric.legendre.helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy_legendre \
                - helmholtz_free_energy_legendre_out \
                + parameters.boltzmann_constant*temperature*np.log(
                    8*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )
            residual_rel = residual_abs/helmholtz_free_energy_legendre
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_helmholtz_free_energy_per_link(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
            end_to_end_length_per_link = \
                model.isotensional.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            helmholtz_free_energy_per_link_legendre = \
                model.isotensional.gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length_per_link
            helmholtz_free_energy_per_link_legendre_out = \
                model.isometric.legendre.helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy_per_link_legendre \
                - helmholtz_free_energy_per_link_legendre_out \
                + parameters.boltzmann_constant*temperature*np.log(
                    8*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )/number_of_links
            residual_rel = residual_abs/helmholtz_free_energy_per_link_legendre
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
            relative_helmholtz_free_energy_legendre = \
                model.isotensional.relative_gibbs_free_energy(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length
            relative_helmholtz_free_energy_legendre_out = \
                model.isometric.legendre.relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_helmholtz_free_energy_legendre \
                - relative_helmholtz_free_energy_legendre_out
            residual_rel = residual_abs/relative_helmholtz_free_energy_legendre
            self.assertLessEqual(
                np.abs(residual_rel),
                3e1 * parameters.rel_tol
            )

    def test_relative_helmholtz_free_energy_per_link(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
            end_to_end_length_per_link = \
                model.isotensional.end_to_end_length_per_link(
                    np.array(force),
                    temperature
                )
            relative_helmholtz_free_energy_per_link_legendre = \
                model.isotensional.relative_gibbs_free_energy_per_link(
                    np.array(force),
                    temperature
                ) + force*end_to_end_length_per_link
            relative_helmholtz_free_energy_per_link_legendre_out = \
                model.isometric.legendre. \
                relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_helmholtz_free_energy_per_link_legendre \
                - relative_helmholtz_free_energy_per_link_legendre_out
            residual_rel = residual_abs / \
                relative_helmholtz_free_energy_per_link_legendre
            self.assertLessEqual(
                np.abs(residual_rel),
                3e1 * parameters.rel_tol
            )

    def test_nondimensional_helmholtz_free_energy(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
            nondimensional_end_to_end_length = \
                model.isotensional.nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                )
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                )
            nondimensional_helmholtz_free_energy_legendre = \
                model.isotensional.nondimensional_gibbs_free_energy(
                    np.array(nondimensional_force),
                    temperature
                ) + nondimensional_force*nondimensional_end_to_end_length
            nondimensional_helmholtz_free_energy_legendre_out = \
                model.isometric.legendre.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy_legendre \
                - nondimensional_helmholtz_free_energy_legendre_out \
                + np.log(
                    8*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )
            residual_rel = residual_abs / \
                nondimensional_helmholtz_free_energy_legendre
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_helmholtz_free_energy_per_link(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
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
                    np.array(nondimensional_force),
                )
            nondimensional_helmholtz_free_energy_per_link_legendre = \
                model.isotensional.nondimensional_gibbs_free_energy_per_link(
                    np.array(nondimensional_force),
                    temperature
                ) + nondimensional_force * \
                nondimensional_end_to_end_length_per_link
            nondimensional_helmholtz_free_energy_per_link_legendre_out = \
                model.isometric.legendre. \
                nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy_per_link_legendre \
                - nondimensional_helmholtz_free_energy_per_link_legendre_out \
                + np.log(
                    8*np.pi**2*hinge_mass*link_length**2 *
                    parameters.boltzmann_constant*temperature /
                    parameters.planck_constant**2
                )/number_of_links
            residual_rel = residual_abs / \
                nondimensional_helmholtz_free_energy_per_link_legendre
            self.assertLessEqual(
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            nondimensional_end_to_end_length = \
                model.isotensional.nondimensional_end_to_end_length(
                    np.array(nondimensional_force),
                )
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                )
            nondimensional_relative_helmholtz_free_energy_legendre = \
                model.isotensional.nondimensional_relative_gibbs_free_energy(
                    np.array(nondimensional_force)
                ) + nondimensional_force*nondimensional_end_to_end_length
            nondimensional_relative_helmholtz_free_energy_legendre_out = \
                model.isometric.legendre. \
                nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_relative_helmholtz_free_energy_legendre \
                - nondimensional_relative_helmholtz_free_energy_legendre_out
            residual_rel = residual_abs / \
                nondimensional_relative_helmholtz_free_energy_legendre
            self.assertLessEqual(
                np.abs(residual_rel),
                3e1 * parameters.rel_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy_per_link(self):
        """Function to test the Legendre transformation
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
            well_width = \
                parameters.well_width_reference + \
                parameters.well_width_scale*(0.5 - np.random.rand())
            model = SWFJC(
                number_of_links,
                link_length,
                hinge_mass,
                well_width
            )
            nondimensional_force = \
                parameters. \
                nondimensional_force_reference + \
                parameters. \
                nondimensional_force_scale * \
                (0.5 - np.random.rand())
            nondimensional_end_to_end_length_per_link = \
                model.isotensional.nondimensional_end_to_end_length_per_link(
                    np.array(nondimensional_force),
                )
            nondim_relative_helmholtz_free_energy_per_link_legendre = \
                model.isotensional. \
                nondimensional_relative_gibbs_free_energy_per_link(
                    np.array(nondimensional_force)
                ) + nondimensional_force * \
                nondimensional_end_to_end_length_per_link
            nondim_relative_helmholtz_free_energy_per_link_legendre_out = \
                model.isometric.legendre. \
                nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondim_relative_helmholtz_free_energy_per_link_legendre \
                - nondim_relative_helmholtz_free_energy_per_link_legendre_out
            residual_rel = residual_abs / \
                nondim_relative_helmholtz_free_energy_per_link_legendre
            self.assertLessEqual(
                np.abs(residual_rel),
                3e1 * parameters.rel_tol
            )
