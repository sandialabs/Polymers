"""Module for dummy tests.

"""

import unittest

from .. import *


class Dummy(unittest.TestCase):
    """Class for dummy tests.

    """
    def test_zero_inverse(self):
        """Function for dummy test.

        """
        self.assertEqual(8, 8)
