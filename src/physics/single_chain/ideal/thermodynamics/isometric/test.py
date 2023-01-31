"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters
from ....test import integrate

parameters = Parameters()
Ideal = physics.single_chain.ideal.thermodynamics.isometric.Ideal


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = Ideal(
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
                Ideal(
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
                Ideal(
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
                Ideal(
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
            model = Ideal(
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


class Normalization(unittest.TestCase):
    """Class for normalization tests.

    """
    def test_equilibrium_distribution(self):
        """Function to test the normalization
        of the equilibrium distribution.

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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )

            def integrand(end_to_end_length):
                return 4.0*np.pi*end_to_end_length**2 * \
                    model.equilibrium_distribution(
                        end_to_end_length
                    )

            integral = integrate(
                integrand,
                parameters.zero*number_of_links*link_length,
                1e1*number_of_links*link_length,
                parameters.points
            )
            self.assertLessEqual(
                np.abs(integral - 1.0),
                parameters.rel_tol
            )

    def test_nondimensional_equilibrium_distribution(self):
        """Function to test the normalization
        of the nondimensional equilibrium distribution.

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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )

            def integrand(nondimensional_end_to_end_length_per_link_per_link):
                return 4.0*np.pi * \
                    nondimensional_end_to_end_length_per_link_per_link**2 * \
                    model.nondimensional_equilibrium_distribution(
                        nondimensional_end_to_end_length_per_link_per_link
                    )

            integral = integrate(
                integrand,
                parameters.zero,
                1e1,
                parameters.points
            )
            self.assertLessEqual(
                np.abs(integral - 1.0),
                parameters.rel_tol
            )

    def test_equilibrium_radial_distribution(self):
        """Function to test the normalization
        of the equilibrium radial distribution.

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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )

            def integrand(end_to_end_length):
                return model.equilibrium_radial_distribution(
                        end_to_end_length
                    )

            integral = integrate(
                integrand,
                parameters.zero*number_of_links*link_length,
                1e1*number_of_links*link_length,
                parameters.points
            )
            self.assertLessEqual(
                np.abs(integral - 1.0),
                parameters.rel_tol
            )

    def test_nondimensional_equilibrium_radial_distribution(self):
        """Function to test the normalization
        of the nondimensional equilibrium radial distribution.

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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )

            def integrand(nondimensional_end_to_end_length_per_link_per_link):
                return model.nondimensional_equilibrium_radial_distribution(
                        nondimensional_end_to_end_length_per_link_per_link
                    )

            integral = integrate(
                integrand,
                parameters.zero,
                1e1,
                parameters.points
            )
            self.assertLessEqual(
                np.abs(integral - 1.0),
                parameters.rel_tol
            )


class Nondimensional(unittest.TestCase):
    """Class for nondimensionalization tests.

    """
    def test_force(self):
        """Function to test the nondimensionalization
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
            model = Ideal(
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
                model.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            force = \
                model.force(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                force / \
                parameters.boltzmann_constant/temperature*link_length \
                - nondimensional_force
            residual_rel = \
                residual_abs / \
                nondimensional_force
            self.assertLessEqual(
                residual_abs,
                parameters.abs_tol
            )
            self.assertLessEqual(
                residual_rel,
                parameters.rel_tol
            )

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
            model = Ideal(
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
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            helmholtz_free_energy = \
                model.helmholtz_free_energy(
                    np.array(end_to_end_length),
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
            model = Ideal(
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
            nondimensional_helmholtz_free_energy_per_link = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            helmholtz_free_energy_per_link = \
                model.helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy_per_link / \
                parameters.boltzmann_constant/temperature \
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
            model = Ideal(
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
            nondimensional_relative_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            relative_helmholtz_free_energy = \
                model.relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_helmholtz_free_energy / \
                parameters.boltzmann_constant/temperature \
                - nondimensional_relative_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_relative_helmholtz_free_energy
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
            model = Ideal(
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
            nondimensional_relative_helmholtz_free_energy_per_link = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            end_to_end_length = nondimensional_end_to_end_length_per_link * \
                number_of_links*link_length
            relative_helmholtz_free_energy_per_link = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                relative_helmholtz_free_energy_per_link / \
                parameters.boltzmann_constant/temperature \
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
            model = Ideal(
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
            helmholtz_free_energy = \
                model.helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_per_link = \
                model.helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
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
            model = Ideal(
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
            relative_helmholtz_free_energy = \
                model.relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy_per_link = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
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
            model = Ideal(
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
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_helmholtz_free_energy_per_link = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
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
            model = Ideal(
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
            nondimensional_relative_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy_per_link = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
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
            model = Ideal(
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
            helmholtz_free_energy = \
                model.helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_0 = \
                model.helmholtz_free_energy(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            relative_helmholtz_free_energy = \
                model.relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy \
                - helmholtz_free_energy_0 \
                - relative_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                relative_helmholtz_free_energy
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
            model = Ideal(
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
            helmholtz_free_energy_per_link = \
                model.helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            helmholtz_free_energy_per_link_0 = \
                model.helmholtz_free_energy_per_link(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            relative_helmholtz_free_energy_per_link = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(end_to_end_length),
                    temperature
                )
            residual_abs = \
                helmholtz_free_energy_per_link \
                - helmholtz_free_energy_per_link_0 \
                - relative_helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                relative_helmholtz_free_energy_per_link
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
            model = Ideal(
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
            nondimensional_helmholtz_free_energy = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_helmholtz_free_energy_0 = \
                model.nondimensional_helmholtz_free_energy(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy \
                - nondimensional_helmholtz_free_energy_0 \
                - nondimensional_relative_helmholtz_free_energy
            residual_rel = \
                residual_abs / \
                nondimensional_relative_helmholtz_free_energy
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
            model = Ideal(
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
            nondimensional_helmholtz_free_energy_per_link = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link),
                    temperature
                )
            nondimensional_helmholtz_free_energy_per_link_0 = \
                model.nondimensional_helmholtz_free_energy_per_link(
                    np.array(parameters.zero),
                    temperature
                )
            nondimensional_relative_helmholtz_free_energy_per_link = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            residual_abs = \
                nondimensional_helmholtz_free_energy_per_link \
                - nondimensional_helmholtz_free_energy_per_link_0 \
                - nondimensional_relative_helmholtz_free_energy_per_link
            residual_rel = \
                residual_abs / \
                nondimensional_relative_helmholtz_free_energy_per_link
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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_helmholtz_free_energy_0 = \
                model.relative_helmholtz_free_energy(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            self.assertLessEqual(
                np.abs(relative_helmholtz_free_energy_0),
                parameters.boltzmann_constant*temperature *
                number_of_links*parameters.zero
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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            temperature = \
                parameters.temperature_reference + \
                parameters.temperature_scale*(0.5 - np.random.rand())
            relative_helmholtz_free_energy_per_link_0 = \
                model.relative_helmholtz_free_energy_per_link(
                    np.array(parameters.zero*number_of_links*link_length),
                    temperature
                )
            self.assertLessEqual(
                np.abs(relative_helmholtz_free_energy_per_link_0),
                parameters.boltzmann_constant*temperature *
                parameters.zero
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the zero
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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_relative_helmholtz_free_energy_0 = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(parameters.zero*number_of_links*link_length)
                )
            self.assertLessEqual(
                np.abs(nondimensional_relative_helmholtz_free_energy_0),
                number_of_links*parameters.zero
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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_relative_helmholtz_free_energy_per_link_0 = \
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(parameters.zero*number_of_links*link_length)
                )
            self.assertLessEqual(
                np.abs(
                    nondimensional_relative_helmholtz_free_energy_per_link_0
                ), parameters.zero
            )

    def test_equilibrium_radial_distribution(self):
        """Function to test the zero
        of the equilibrium radial distribution.

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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            equilibrium_radial_distribution_0 = \
                model.equilibrium_radial_distribution(
                    np.array(parameters.zero*number_of_links*link_length)
                )
            self.assertLessEqual(
                np.abs(equilibrium_radial_distribution_0),
                parameters.zero
            )

    def test_nondimensional_equilibrium_radial_distribution(self):
        """Function to test the zero
        of the nondimensional equilibrium radial distribution.

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
            model = Ideal(
                number_of_links,
                link_length,
                hinge_mass
            )
            nondimensional_equilibrium_radial_distribution_0 = \
                model.nondimensional_equilibrium_radial_distribution(
                    np.array(parameters.zero)
                )
            self.assertLessEqual(
                np.abs(nondimensional_equilibrium_radial_distribution_0),
                parameters.zero
            )


class Connection(unittest.TestCase):
    """Class for connection tests.

    """
    def test_force(self):
        """Function to test the connection
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
            model = Ideal(
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
                model.force(
                    np.array(end_to_end_length),
                    temperature
                )
            h_step = parameters.rel_tol * \
                number_of_links*link_length
            force_from_derivative = (
                model.relative_helmholtz_free_energy(
                    np.array(end_to_end_length + 0.5*h_step),
                    temperature
                )
                - model.relative_helmholtz_free_energy(
                    np.array(end_to_end_length - 0.5*h_step),
                    temperature
                ))/h_step
            residual_abs = \
                force \
                - force_from_derivative
            residual_rel = residual_abs/force
            self.assertLessEqual(
                np.abs(residual_rel), h_step
            )

    def test_nondimensional_force(self):
        """Function to test the connection
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
            model = Ideal(
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
                model.nondimensional_force(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            h_step = parameters.rel_tol * \
                number_of_links*link_length
            nondimensional_force_from_derivative = (
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(
                        nondimensional_end_to_end_length_per_link + 0.5*h_step
                    )
                )
                - model.nondimensional_relative_helmholtz_free_energy_per_link(
                    np.array(
                        nondimensional_end_to_end_length_per_link - 0.5*h_step
                    )
                ))/h_step
            residual_abs = \
                nondimensional_force \
                - nondimensional_force_from_derivative
            residual_rel = residual_abs/nondimensional_force
            self.assertLessEqual(
                np.abs(residual_rel), h_step
            )

    def test_relative_helmholtz_free_energy(self):
        """Function to test the connection
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
            model = Ideal(
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
            relative_helmholtz_free_energy = \
                model.relative_helmholtz_free_energy(
                    np.array(end_to_end_length),
                    temperature
                )
            relative_helmholtz_free_energy_from_connection = \
                parameters.boltzmann_constant*temperature * np.log(
                    model.equilibrium_distribution(np.array(
                        parameters.zero*number_of_links*link_length
                    )) /
                    model.equilibrium_distribution(np.array(
                        end_to_end_length
                    ))
                )
            residual_abs = \
                relative_helmholtz_free_energy \
                - relative_helmholtz_free_energy_from_connection
            residual_rel = residual_abs/relative_helmholtz_free_energy
            self.assertLessEqual(
                np.abs(residual_rel), parameters.rel_tol
            )

    def test_nondimensional_relative_helmholtz_free_energy(self):
        """Function to test the connection
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
            model = Ideal(
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
            nondimensional_relative_helmholtz_free_energy = \
                model.nondimensional_relative_helmholtz_free_energy(
                    np.array(nondimensional_end_to_end_length_per_link)
                )
            nondimensional_relative_helmholtz_free_energy_from_connection = \
                np.log(
                    model.nondimensional_equilibrium_distribution(np.array(
                        parameters.zero
                    )) /
                    model.nondimensional_equilibrium_distribution(np.array(
                        nondimensional_end_to_end_length_per_link
                    ))
                )
            residual_abs = \
                nondimensional_relative_helmholtz_free_energy \
                - nondimensional_relative_helmholtz_free_energy_from_connection
            residual_rel = residual_abs / \
                nondimensional_relative_helmholtz_free_energy
            self.assertLessEqual(
                np.abs(residual_rel), parameters.rel_tol
            )
