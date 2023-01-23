"""Module to test ...

"""

import unittest
from polymers import physics
from ..test import Parameters

parameters = Parameters()
FJC = physics.single_chain.fjc.thermodynamics.isotensional.legendre.FJC


class Base(unittest.TestCase):
    """Class to test ...

    """
    def test_init(self):
        """Function to test ...

        """
        _ = FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference
        )
