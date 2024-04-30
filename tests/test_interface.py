import os.path
import unittest
import sys

# add parent directory to sys.path
current_script_path = os.path.abspath(__file__)  # Absolute path to the current script
current_dir = os.path.dirname(current_script_path)  # Directory containing the current script
parent_dir = os.path.dirname(current_dir)  # Parent directory
src_dir = os.path.join(parent_dir, "src") # Source file directory
sys.path.insert(0, parent_dir)

import plot_orbitals


class TestPyCppInterface(unittest.TestCase):
    def test_generate_points(self):
        H2_points = plot_orbitals.generate_h2_points()
        self.assertTrue(len(H2_points) > 40)