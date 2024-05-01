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
        self.assertTrue(os.path.exists("H0._0.000_0.000_0.000.txt"))
        self.assertTrue(os.path.exists("H1s_1.398_0.000_0.000.txt"))


    def test_read_orbital_points(self):
        target_num_points = 40 
        # 40 super arbitrary
        # seemed like OK starting point balanced density/sparsity for interpolation
        self.assertTrue(len(H2_points) > target_num_points) 