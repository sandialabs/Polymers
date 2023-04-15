"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters
from ...test import integrate

parameters = Parameters()
FJC = physics.single_chain.fjc.thermodynamics.FJC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = FJC(
                parameters.number_of_links_minimum,
                parameters.link_length_reference,
                parameters.hinge_mass_reference
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
                FJC(
                    number_of_links,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference
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
                FJC(
                    parameters.number_of_links_minimum,
                    link_length,
                    parameters.hinge_mass_reference
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
                FJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass
                ).hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                np.abs(residual_abs),
                parameters.abs_tol
            )
            self.assertLessEqual(
                np.abs(residual_rel),
                parameters.rel_tol
            )


class ThermodynamicLimit(unittest.TestCase):
    """Class for thermodynamic limit tests.

    """
    def test_end_to_end_length(self):
        """Function to test the thermodynamic limit
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
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - helmholtz_free_energy_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - helmholtz_free_energy_per_link_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - nondimensional_helmholtz_free_energy_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - nondimensional_helmholtz_free_energy_per_link_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - gibbs_free_energy_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - gibbs_free_energy_per_link_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - nondimensional_gibbs_free_energy_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
                - nondimensional_gibbs_free_energy_per_link_out
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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
            number_of_links = parameters.number_of_links_maximum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
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


class ModifiedCanonicalStrongPotentialIsometric(unittest.TestCase):
    """Class for the modified canonical ensemble
    becoming the isometric ensemble under strong potential tests.

    """
    def test_force(self):
        """Function to test the behavior
        of the force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.
                        force(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        force(
                            end_to_end_length,
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.modified_canonical.
                        force(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_force(self):
        """Function to test the behavior
        of the nondimensional force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.
                        nondimensional_force(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        ) -
                        model.isometric.
                        nondimensional_force(
                            nondimensional_end_to_end_length_per_link
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.
                        nondimensional_force(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the behavior
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.
                        relative_helmholtz_free_energy(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        relative_helmholtz_free_energy(
                            end_to_end_length,
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.modified_canonical.
                        relative_helmholtz_free_energy(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_helmholtz_free_energy_per_link(self):
        """Function to test the behavior
        of the relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.
                        relative_helmholtz_free_energy_per_link(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        relative_helmholtz_free_energy_per_link(
                            end_to_end_length,
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.modified_canonical.
                        relative_helmholtz_free_energy_per_link(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the behavior
        of the nondimensional relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.
                        nondimensional_relative_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        ) -
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.
                        nondimensional_relative_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy_per_link(self):
        """Function to test the behavior
        of the nondimensional relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        ) -
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )


class ModifiedCanonicalWeakPotentialIsotensional(unittest.TestCase):
    """Class for the modified canonical ensemble
    becoming the isotensional ensemble under weak potential tests.

    """
    def test_end_to_end_length(self):
        """Function to test the behavior
        of the end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical.force(
                        potential_distance,
                        potential_stiffness,
                        temperature
                    )
                    return (
                        model.modified_canonical.
                        end_to_end_length(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        end_to_end_length(
                            force,
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.
                        end_to_end_length(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the behavior
        of the end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical.force(
                        potential_distance,
                        potential_stiffness,
                        temperature
                    )
                    return (
                        model.modified_canonical.
                        end_to_end_length_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        end_to_end_length_per_link(
                            force,
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.
                        end_to_end_length_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the behavior
        of the nondimensional end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.
                        nondimensional_end_to_end_length(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_end_to_end_length(
                            nondimensional_force
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.
                        nondimensional_end_to_end_length(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the behavior
        of the nondimensional end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_force
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_gibbs_free_energy(self):
        """Function to test the behavior
        of the Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical.force(
                    np.array(potential_distance_ref),
                    potential_stiffness,
                    temperature
                )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical.force(
                        potential_distance,
                        potential_stiffness,
                        temperature
                    )
                    return (
                        model.modified_canonical.
                        gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        gibbs_free_energy(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        gibbs_free_energy(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.
                        gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical.force(
                    np.array(potential_distance_ref),
                    potential_stiffness,
                    temperature
                )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical.force(
                        potential_distance,
                        potential_stiffness,
                        temperature
                    )
                    return (
                        model.modified_canonical.
                        gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        gibbs_free_energy_per_link(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        gibbs_free_energy_per_link(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.
                        gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_gibbs_free_energy(self):
        """Function to test the behavior
        of the relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical.force(
                    np.array(potential_distance_ref),
                    potential_stiffness,
                    temperature
                )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical.force(
                        potential_distance,
                        potential_stiffness,
                        temperature
                    )
                    return (
                        model.modified_canonical.
                        relative_gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        relative_gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        relative_gibbs_free_energy(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        relative_gibbs_free_energy(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.
                        relative_gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        relative_gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical.force(
                    np.array(potential_distance_ref),
                    potential_stiffness,
                    temperature
                )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical.force(
                        potential_distance,
                        potential_stiffness,
                        temperature
                    )
                    return (
                        model.modified_canonical.
                        relative_gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        relative_gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        relative_gibbs_free_energy_per_link(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        relative_gibbs_free_energy_per_link(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.
                        relative_gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        relative_gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_gibbs_free_energy(self):
        """Function to test the behavior
        of the nondimensional Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        nondimensional_gibbs_free_energy(
                            nondimensional_force,
                            temperature
                        ) +
                        model.isotensional.
                        nondimensional_gibbs_free_energy(
                            np.array(nondimensional_force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the nondimensional Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        nondimensional_gibbs_free_energy_per_link(
                            nondimensional_force,
                            temperature
                        ) +
                        model.isotensional.
                        nondimensional_gibbs_free_energy_per_link(
                            np.array(nondimensional_force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.
                        nondimensional_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the behavior
        of the nondimensional relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy(
                            nondimensional_force
                        ) +
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy(
                            np.array(nondimensional_force_ref)
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the nondimensional relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            nondimensional_force
                        ) +
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            np.array(nondimensional_force_ref)
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )


class ModifiedCanonicalAsymptoticStrongPotentialIsometric(unittest.TestCase):
    """Class for the modified canonical ensemble asymptotics
    becoming the isometric ensemble under strong potential tests.

    """
    def test_force(self):
        """Function to test the behavior
        of the force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        force(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        force(
                            end_to_end_length,
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.isometric.
                        force(
                            end_to_end_length,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_force(self):
        """Function to test the behavior
        of the nondimensional force.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_force(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        ) -
                        model.isometric.
                        nondimensional_force(
                            nondimensional_end_to_end_length_per_link
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.isometric.
                        nondimensional_force(
                            nondimensional_end_to_end_length_per_link
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_helmholtz_free_energy(self):
        """Function to test the behavior
        of the Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        helmholtz_free_energy(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        helmholtz_free_energy(
                            end_to_end_length,
                            temperature
                        ) +
                        model.isometric.
                        helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.isometric.
                        helmholtz_free_energy(
                            end_to_end_length,
                            temperature
                        ) -
                        model.isometric.
                        helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_helmholtz_free_energy_per_link(self):
        """Function to test the behavior
        of the Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        helmholtz_free_energy_per_link(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        helmholtz_free_energy_per_link(
                            end_to_end_length,
                            temperature
                        ) +
                        model.isometric.
                        helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.isometric.
                        helmholtz_free_energy_per_link(
                            end_to_end_length,
                            temperature
                        ) -
                        model.isometric.
                        helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the behavior
        of the relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        relative_helmholtz_free_energy(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        relative_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        relative_helmholtz_free_energy(
                            end_to_end_length,
                            temperature
                        ) +
                        model.isometric.
                        relative_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.isometric.
                        relative_helmholtz_free_energy(
                            end_to_end_length,
                            temperature
                        ) -
                        model.isometric.
                        relative_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_helmholtz_free_energy_per_link(self):
        """Function to test the behavior
        of the relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(end_to_end_length):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        relative_helmholtz_free_energy_per_link(
                            end_to_end_length,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        relative_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        relative_helmholtz_free_energy_per_link(
                            end_to_end_length,
                            temperature
                        ) +
                        model.isometric.
                        relative_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                def integrand_denominator(end_to_end_length):
                    return (
                        model.isometric.
                        relative_helmholtz_free_energy_per_link(
                            end_to_end_length,
                            temperature
                        ) -
                        model.isometric.
                        relative_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small *
                                number_of_links*link_length
                            ),
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero*number_of_links*link_length,
                    parameters.nondimensional_potential_distance_small *
                    number_of_links*link_length,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_helmholtz_free_energy(self):
        """Function to test the behavior
        of the nondimensional Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        nondimensional_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link,
                            temperature
                        ) +
                        model.isometric.
                        nondimensional_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            temperature
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.isometric.
                        nondimensional_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link,
                            temperature
                        ) -
                        model.isometric.
                        nondimensional_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_helmholtz_free_energy_per_link(self):
        """Function to test the behavior
        of the nondimensional Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.isometric.
                        nondimensional_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link,
                            temperature
                        ) +
                        model.isometric.
                        nondimensional_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            temperature
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.isometric.
                        nondimensional_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link,
                            temperature
                        ) -
                        model.isometric.
                        nondimensional_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the behavior
        of the nondimensional relative Helmholtz free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_relative_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_relative_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            nondimensional_potential_stiffness
                        ) -
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link
                        ) +
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            )
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy(
                            nondimensional_end_to_end_length_per_link
                        ) -
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            )
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy_per_link(self):
        """Function to test the behavior
        of the nondimensional relative Helmholtz free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.asymptotic.strong_potential.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            ),
                            nondimensional_potential_stiffness
                        ) -
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link
                        ) +
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            )
                        )
                    )**2

                def integrand_denominator(
                    nondimensional_end_to_end_length_per_link
                ):
                    return (
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            nondimensional_end_to_end_length_per_link
                        ) -
                        model.isometric.
                        nondimensional_relative_helmholtz_free_energy_per_link(
                            np.array(
                                parameters.
                                nondimensional_potential_distance_small
                            )
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.zero,
                    parameters.nondimensional_potential_distance_small,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_large
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_large *
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )


class ModifiedCanonicalAsymptoticWeakPotentialIsotensional(unittest.TestCase):
    """Class for the modified canonical ensemble asymptotics
    becoming the isotensional ensemble under weak potential tests.

    """
    def test_end_to_end_length(self):
        """Function to test the behavior
        of the end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical. \
                        asymptotic.weak_potential. \
                        force(
                            potential_distance,
                            potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        end_to_end_length(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        end_to_end_length(
                            force,
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        end_to_end_length(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_end_to_end_length_per_link(self):
        """Function to test the behavior
        of the end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical. \
                        asymptotic.weak_potential. \
                        force(
                            potential_distance,
                            potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        end_to_end_length_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        end_to_end_length_per_link(
                            force,
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        end_to_end_length_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length(self):
        """Function to test the behavior
        of the nondimensional end-to-end length.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.asymptotic.weak_potential. \
                        nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_end_to_end_length(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_end_to_end_length(
                            nondimensional_force
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_end_to_end_length(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_end_to_end_length_per_link(self):
        """Function to test the behavior
        of the nondimensional end-to-end length per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )

            def residual_rel(nondimensional_potential_stiffness):

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.asymptotic.weak_potential. \
                        nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_force
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_end_to_end_length_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_gibbs_free_energy(self):
        """Function to test the behavior
        of the Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical. \
                    asymptotic.weak_potential. \
                    force(
                        np.array(potential_distance_ref),
                        potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical. \
                        asymptotic.weak_potential. \
                        force(
                            potential_distance,
                            potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        gibbs_free_energy(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        gibbs_free_energy(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical. \
                    asymptotic.weak_potential. \
                    force(
                        np.array(potential_distance_ref),
                        potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical. \
                        asymptotic.weak_potential. \
                        force(
                            potential_distance,
                            potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        gibbs_free_energy_per_link(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        gibbs_free_energy_per_link(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_gibbs_free_energy(self):
        """Function to test the behavior
        of the relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical. \
                    asymptotic.weak_potential. \
                    force(
                        np.array(potential_distance_ref),
                        potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical. \
                        asymptotic.weak_potential. \
                        force(
                            potential_distance,
                            potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        relative_gibbs_free_energy(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        relative_gibbs_free_energy(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_relative_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )*number_of_links*link_length

            def residual_rel(nondimensional_potential_stiffness):
                potential_stiffness = nondimensional_potential_stiffness / \
                    link_length**2 * \
                    parameters.boltzmann_constant*temperature
                force_ref = model.modified_canonical. \
                    asymptotic.weak_potential. \
                    force(
                        np.array(potential_distance_ref),
                        potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    force = model.modified_canonical. \
                        asymptotic.weak_potential. \
                        force(
                            potential_distance,
                            potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        relative_gibbs_free_energy_per_link(
                            force,
                            temperature
                        ) +
                        model.isotensional.
                        relative_gibbs_free_energy_per_link(
                            np.array(force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    potential_distance = nondimensional_potential_distance * \
                        number_of_links*link_length
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy_per_link(
                            potential_distance,
                            potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        relative_gibbs_free_energy_per_link(
                            np.array(potential_distance_ref),
                            potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_gibbs_free_energy(self):
        """Function to test the behavior
        of the nondimensional Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.asymptotic.weak_potential. \
                    nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.asymptotic.weak_potential. \
                        nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        nondimensional_gibbs_free_energy(
                            nondimensional_force,
                            temperature
                        ) +
                        model.isotensional.
                        nondimensional_gibbs_free_energy(
                            np.array(nondimensional_force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the nondimensional Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.asymptotic.weak_potential. \
                    nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.asymptotic.weak_potential. \
                        nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.isotensional.
                        nondimensional_gibbs_free_energy_per_link(
                            nondimensional_force,
                            temperature
                        ) +
                        model.isotensional.
                        nondimensional_gibbs_free_energy_per_link(
                            np.array(nondimensional_force_ref),
                            temperature
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness,
                            temperature
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness,
                            temperature
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_gibbs_free_energy(self):
        """Function to test the behavior
        of the nondimensional relative Gibbs free energy.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.asymptotic.weak_potential. \
                    nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.asymptotic.weak_potential. \
                        nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy(
                            nondimensional_force
                        ) +
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy(
                            np.array(nondimensional_force_ref)
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )

    def test_nondimensional_relative_gibbs_free_energy_per_link(self):
        """Function to test the behavior
        of the nondimensional relative Gibbs free energy per link.

        """
        for _ in range(parameters.number_of_loops):
            number_of_links = \
                parameters.number_of_links_maximum - \
                parameters.number_of_links_minimum
            link_length = \
                parameters.link_length_reference + \
                parameters.link_length_scale*(0.5 - np.random.rand())
            hinge_mass = \
                parameters.hinge_mass_reference + \
                parameters.hinge_mass_scale*(0.5 - np.random.rand())
            model = FJC(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_potential_distance_ref = \
                np.array(
                    parameters.nondimensional_potential_distance_large_1
                )

            def residual_rel(nondimensional_potential_stiffness):
                nondimensional_force_ref = \
                    model.modified_canonical.asymptotic.weak_potential. \
                    nondimensional_force(
                        np.array(nondimensional_potential_distance_ref),
                        nondimensional_potential_stiffness
                    )

                def integrand_numerator(nondimensional_potential_distance):
                    nondimensional_force = \
                        model.modified_canonical.asymptotic.weak_potential. \
                        nondimensional_force(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        )
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        ) -
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            nondimensional_force
                        ) +
                        model.isotensional.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            np.array(nondimensional_force_ref)
                        )
                    )**2

                def integrand_denominator(nondimensional_potential_distance):
                    return (
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            nondimensional_potential_distance,
                            nondimensional_potential_stiffness
                        ) -
                        model.modified_canonical.asymptotic.weak_potential.
                        nondimensional_relative_gibbs_free_energy_per_link(
                            np.array(nondimensional_potential_distance_ref),
                            nondimensional_potential_stiffness
                        )
                    )**2

                numerator = integrate(
                    integrand_numerator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                denominator = integrate(
                    integrand_denominator,
                    parameters.nondimensional_potential_distance_large_1,
                    parameters.nondimensional_potential_distance_large_2,
                    parameters.points)
                return np.sqrt(numerator/denominator)

            residual_rel_1 = residual_rel(
                parameters.nondimensional_potential_stiffness_small
            )
            residual_rel_2 = residual_rel(
                parameters.nondimensional_potential_stiffness_small /
                parameters.log_log_scale
            )
            log_log_slope = np.log(residual_rel_2/residual_rel_1) / \
                np.log(parameters.log_log_scale)
            self.assertLessEqual(
                np.abs(log_log_slope + 1.0),
                parameters.log_log_tol
            )
