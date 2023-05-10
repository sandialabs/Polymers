"""Module to test the local module.

"""
import unittest
import numpy as np
from polymers import physics
from ..test import Parameters

parameters = Parameters()
CUFJC = physics.single_chain.ufjc.composite. \
    thermodynamics.isometric.asymptotic.reduced.legendre.CUFJC


class Base(unittest.TestCase):
    """Class for basic tests.

    """
    def test_init(self):
        """Function to test instantiation.

        """
        for _ in range(parameters.number_of_loops):
            _ = CUFJC(
                parameters.number_of_links_minimum,
                parameters.link_length_reference,
                parameters.hinge_mass_reference,
                parameters.number_of_bonds_minimum,
                parameters.bond_stiffness_reference,
                parameters.bond_energy_reference,
                parameters.bond_scission_energy_reference,
                parameters.bond_attempt_frequency_reference
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
                CUFJC(
                    number_of_links,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.number_of_bonds_minimum,
                    parameters.bond_stiffness_reference,
                    parameters.bond_energy_reference,
                    parameters.bond_scission_energy_reference,
                    parameters.bond_attempt_frequency_reference
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
                CUFJC(
                    parameters.number_of_links_minimum,
                    link_length,
                    parameters.hinge_mass_reference,
                    parameters.number_of_bonds_minimum,
                    parameters.bond_stiffness_reference,
                    parameters.bond_energy_reference,
                    parameters.bond_scission_energy_reference,
                    parameters.bond_attempt_frequency_reference
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
                CUFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass,
                    parameters.number_of_bonds_minimum,
                    parameters.bond_stiffness_reference,
                    parameters.bond_energy_reference,
                    parameters.bond_scission_energy_reference,
                    parameters.bond_attempt_frequency_reference
                ).hinge_mass
            )

    def test_number_of_bonds(self):
        """Function to test the number of bonds during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            number_of_bonds = \
                np.random.randint(
                    parameters.number_of_bonds_minimum,
                    high=parameters.number_of_bonds_maximum
                )
            self.assertEqual(
                number_of_bonds,
                CUFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    number_of_bonds,
                    parameters.bond_stiffness_reference,
                    parameters.bond_energy_reference,
                    parameters.bond_scission_energy_reference,
                    parameters.bond_attempt_frequency_reference
                ).number_of_bonds
            )

    def test_bond_stiffness(self):
        """Function to test the bond stiffness during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            bond_stiffness = \
                parameters.bond_stiffness_reference + \
                parameters.bond_stiffness_scale*(0.5 - np.random.rand())
            self.assertEqual(
                bond_stiffness,
                CUFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.number_of_bonds_minimum,
                    bond_stiffness,
                    parameters.bond_energy_reference,
                    parameters.bond_scission_energy_reference,
                    parameters.bond_attempt_frequency_reference
                ).bond_stiffness
            )

    def test_bond_energy(self):
        """Function to test the bond energy during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            bond_energy = \
                parameters.bond_energy_reference + \
                parameters.bond_energy_scale*(0.5 - np.random.rand())
            self.assertEqual(
                bond_energy,
                CUFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.number_of_bonds_minimum,
                    parameters.bond_stiffness_reference,
                    bond_energy,
                    parameters.bond_scission_energy_reference,
                    parameters.bond_attempt_frequency_reference
                ).bond_energy
            )

    def test_bond_scission_energy(self):
        """Function to test the bond scission energy during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            bond_scission_energy = \
                parameters.bond_scission_energy_reference + \
                parameters.bond_scission_energy_scale*(0.5 - np.random.rand())
            self.assertEqual(
                bond_scission_energy,
                CUFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.number_of_bonds_minimum,
                    parameters.bond_stiffness_reference,
                    parameters.bond_energy_reference,
                    bond_scission_energy,
                    parameters.bond_attempt_frequency_reference
                ).bond_scission_energy
            )

    def test_bond_attempt_frequency(self):
        """Function to test the bond attempt frequency during instantiation.

        """
        for _ in range(parameters.number_of_loops):
            bond_attempt_frequency = \
                parameters.bond_attempt_frequency_reference + \
                parameters.bond_attempt_frequency_scale*(
                    0.5 - np.random.rand()
                )
            self.assertEqual(
                bond_attempt_frequency,
                CUFJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                    parameters.number_of_bonds_minimum,
                    parameters.bond_stiffness_reference,
                    parameters.bond_energy_reference,
                    parameters.bond_scission_energy_reference,
                    bond_attempt_frequency
                ).bond_attempt_frequency
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
            number_of_bonds = \
                np.random.randint(
                    parameters.number_of_bonds_minimum,
                    high=parameters.number_of_bonds_maximum
                )
            bond_stiffness = \
                parameters.bond_stiffness_reference + \
                parameters.bond_stiffness_scale*(0.5 - np.random.rand())
            bond_energy = \
                parameters.bond_energy_reference + \
                parameters.bond_energy_scale*(0.5 - np.random.rand())
            bond_scission_energy = \
                parameters.bond_scission_energy_reference + \
                parameters.bond_scission_energy_scale*(0.5 - np.random.rand())
            bond_attempt_frequency = \
                parameters.bond_attempt_frequency_reference + \
                parameters.bond_attempt_frequency_scale*(
                    0.5 - np.random.rand()
                )
            model = CUFJC(
                number_of_links,
                link_length,
                hinge_mass,
                number_of_bonds,
                bond_stiffness,
                bond_energy,
                bond_scission_energy,
                bond_attempt_frequency
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
                number_of_bonds,
                model.number_of_bonds
            )
            self.assertEqual(
                bond_stiffness,
                model.bond_stiffness
            )
            self.assertEqual(
                bond_energy,
                model.bond_energy
            )
            self.assertEqual(
                bond_scission_energy,
                model.bond_scission_energy
            )
            self.assertEqual(
                bond_attempt_frequency,
                model.bond_attempt_frequency
            )
