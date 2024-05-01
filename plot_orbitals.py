import os.path

import subprocess # subprocess is for calling the c++ calculations, see generate_molecule_points

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata

import current_plot_test


def generate_atomic_orbital_points(filename: str):
    input_path = os.path.join("sample_input", filename)
    assert os.path.exists(input_path), "File not found: May need to run from top-level dir until cleanup"

    assert os.path.exists(os.path.join("cpp", "mo_points_main")), "C++ Executable not compiled"
    command = f"./mo_points_main {input_path}" # mo_points_main takes the path to the file as an arg
    subprocess.run(command, shell=True, check=True)
    
    # Check if the output file was created
    assert os.path.exists(output_path), "Output file not created by the C++ program"
    
    return output_path

def generate_h2_points():
    A = np.loadtxt("cpp_src/s_0.00_0.00_0.00.txt")
    A = A.reshape(50, 50, 50)

    B = np.loadtxt("cpp_src/s_1.40_0.00_0.00.txt")
    B = B.reshape(50, 50, 50)

    return [A, B]


if __name__ == "__main__":
    test_points = generate_h2_points()
    level = np.max(test_points[0]) * 0.5
    colors = ["blue", "red"]

    current_plot_test.plot_filtered_surfaces(test_points, level, colors, [(0.0, 0.0, 0.0), (1.3984, 0.0, 0.0)], (-5.0, 5.0))