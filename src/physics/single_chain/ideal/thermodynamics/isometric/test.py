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
