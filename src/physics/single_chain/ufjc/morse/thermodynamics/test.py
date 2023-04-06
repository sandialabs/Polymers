"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters

parameters = Parameters()
MORSEFJC = physics.single_chain.ufjc. \
    morse.thermodynamics.MORSEFJC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = MORSEFJC(
                parameters.number_of_links_minimum,
                parameters.link_length_reference,
                parameters.hinge_mass_reference,
                parameters.link_stiffness_reference,
                parameters.link_energy_reference
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
                MORSEFJC(
                    number_of_links,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.link_stiffness_reference,
                    parameters.link_energy_reference
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
                MORSEFJC(
                    parameters.number_of_links_minimum,
                    link_length,
                    parameters.hinge_mass_reference,
                    parameters.link_stiffness_reference,
                    parameters.link_energy_reference
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
                MORSEFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass,
                    parameters.link_stiffness_reference,
                    parameters.link_energy_reference
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
                MORSEFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    link_stiffness,
                    parameters.link_energy_reference
                ).link_stiffness
            )

    def test_link_energy(self):
        """Function to test the link energy during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            link_energy = \
                parameters.link_energy_reference + \
                parameters.link_energy_scale*(0.5 - np.random.rand())
            self.assertEqual(
                link_energy,
                MORSEFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.link_stiffness_reference,
                    link_energy,
                ).link_energy
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
            link_energy = \
                parameters.link_energy_reference + \
                parameters.link_energy_scale*(0.5 - np.random.rand())
            model = MORSEFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
                link_energy
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
            self.assertEqual(
                link_energy,
                model.link_energy
            )
