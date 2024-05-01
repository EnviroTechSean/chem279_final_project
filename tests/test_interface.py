import os.path
import unittest
import sys

# add parent directory to sys.path
# usage: pytest -v from this or parent directory
current_script_path = os.path.abspath(__file__)  # Absolute path to the current script
current_dir = os.path.dirname(current_script_path)  # Directory containing the current script
parent_dir = os.path.dirname(current_dir)  # Parent directory
src_dir = os.path.join(parent_dir, "src") # Source file directory
sys.path.insert(0, parent_dir)

import plot_orbitals


class TestPyCppInterface(unittest.TestCase):
    def test_generate_atomic_orbital_points(self):
        plot_orbitals.generate_atomic_orbital_points("H2.txt")
        self.assertTrue(os.path.exists("Hs_0.000_0.000_0.000.txt"))
        self.assertTrue(os.path.exists("Hs_1.398_0.000_0.000.txt"))


    def test_read_orbital_points(self):
        plot_orbitals.generate_atomic_orbital_points("H2.txt")
        target_num_points = 40 # 40 is arbitrary guess at sparse/dense balance for interpolation
        H2_points = plot_orbitals.read_orbital_points("H2.txt")
        self.assertTrue(len(H2_points) > target_num_points) 