"""Module for dummy tests.

"""

import unittest


class Dummy(unittest.TestCase):
    """Class for dummy tests.

    """
    def test_dummy(self):
        """Function for dummy test.

        """
        self.assertEqual(8, 8)
