import os
import os.path

import subprocess # subprocess is for calling the c++ calculations, see generate_molecule_points

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata

import current_plot_development

def generate_atomic_orbital_points(filename: str):
    input_path = os.path.join("sample_input", filename)
    assert os.path.exists(input_path), "File not found: May need to run from top-level dir until cleanup"

    cwd = os.getcwd()
    executable_path = "./mo_points_main"
    assert os.path.exists(executable_path), "C++ Executable not compiled"
    input_path = os.path.join("../", input_path) # Adjust input filename per upcoming cwd change
    command = f"{executable_path} {input_path}"  # mo_points_main takes the path to the file as an arg
    subprocess.run(command, shell=True, check=True, cwd=os.path.join(cwd, "cpp"))
    
    os.chdir(cwd)  # Change back to the original directory

def generate_atomic_orbital_points(filename: str):
    input_path = os.path.join("sample_input", filename)
    assert os.path.exists(input_path), "File not found: May need to run from top-level dir until cleanup"

    executable_path = os.path.join("cpp", "mo_points_main")
    assert os.path.exists(executable_path), "C++ Executable not compiled"
    command = f"{executable_path} {input_path}" # mo_points_main takes the path to the file as an arg
    subprocess.run(command, shell=True, check=True)
    
    return output_path

def read_h2_points():
    A = np.loadtxt("cpp/Hs000_0.000_0.000_0.000.txt")
    A = A.reshape(50, 50, 50)
    B = np.loadtxt("cpp/Hs000_1.398_0.000_0.000.txt")
    B = B.reshape(50, 50, 50)
    return [A, B]


def read_ho_points():
    A = np.loadtxt("cpp/Hs000_0.000_0.000_0.000.txt")
    A = A.reshape(50, 50, 50)
    B = np.loadtxt("cpp/Op100_1.830_0.000_0.000.txt")
    B = B.reshape(50, 50, 50)
    C = np.loadtxt("cpp/Op010_1.830_0.000_0.000.txt")
    C = C.reshape(50, 50, 50)
    D = np.loadtxt("cpp/Op001_1.830_0.000_0.000.txt")
    D = D.reshape(50, 50, 50)
    E = np.loadtxt("cpp/Os000_1.830_0.000_0.000.txt")
    E = E.reshape(50, 50, 50)
    return [A, B, C, D, E]
    # return [A, E]

    # A = np.loadtxt("cpp/Hs000_0.000_0.000_0.000.txt")
    # A = A.reshape(50, 50, 50)
    # B = np.loadtxt("cpp/Op100_1.830_0.000_0.000.txt")
    # B = B.reshape(50, 50, 50)
    # C = np.loadtxt("cpp/Op010_1.830_0.000_0.000.txt")
    # C = C.reshape(50, 50, 50)
    # D = np.loadtxt("cpp/Op001_1.830_0.000_0.000.txt")
    # D = D.reshape(50, 50, 50)
    # E = np.loadtxt("cpp/Os000_1.830_0.000_0.000.txt")
    # E = E.reshape(50, 50, 50)
    # return [B, C, D]




if __name__ == "__main__":
    # test_points = read_h2_points()
    test_points = read_ho_points()

    s_level = np.max(test_points[0]) * 0.5
    p_level = np.max(test_points[1]) * 0.3
    levels = [s_level, p_level, p_level, p_level, s_level]
    colors = ["blue", "red", "green", "yellow", "blue"]
    # colors = ["red", "blue", "green"]
    # colors = ["red", "blue"]

    current_plot_development.plot_filtered_surfaces(test_points, levels, colors, [(0.0, 0.0, 0.0), (1.830, 0.0, 0.0), (1.830, 0.0, 0.0), (1.830, 0.0, 0.0), (1.830, 0.0, 0.0)], (-5.0, 5.0))
    # current_plot_development.plot_filtered_surfaces(test_points, level, colors, [(0.0, 0.0, 0.0), (1.830, 0.0, 0.0)], (-5.0, 5.0))
    # current_plot_development.plot_point_cloud(test_points, (-5.0, 5.0))

    # generate_atomic_orbital_points("H2.txt")