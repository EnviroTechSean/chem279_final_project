import os
import os.path

import subprocess # subprocess is for calling the c++ calculations, see generate_molecule_points

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata

import current_plot_development

def read_orbital_points(filepath):
    assert os.path.exists(filepath)
    data = np.loadtxt(filepath)
    data = data.reshape(50, 50, 50)
    return data

def plot_ho_atomic_orbitals(centers):
    H1s = read_orbital_points("cpp/Hs000_0.000_0.000_0.000.txt")
    O1px = read_orbital_points("cpp/Op100_1.830_0.000_0.000.txt")
    O1py = read_orbital_points("cpp/Op010_1.830_0.000_0.000.txt")
    O1pz = read_orbital_points("cpp/Op001_1.830_0.000_0.000.txt")
    O2s = read_orbital_points("cpp/Os000_1.830_0.000_0.000.txt")

    # Note the redundant p orbitals because of the two different signed isosurfaces
    # This is a hack of course: there is a cleaner factoring that I haven't come up with yet
    test_points = [H1s, O1px, O1py, O1pz, O2s]

    s_level = np.max(test_points[0]) * 0.5
    p_level = np.max(test_points[1]) * 0.4
    levels = [s_level, p_level, p_level, p_level, s_level]

    colors = ["blue", "red", "green", "yellow", "blue"]

    current_plot_development.plot_filtered_surfaces(test_points, levels, colors, centers, (-10.0, 10.0), True)

    # generate_atomic_orbital_points("H2.txt")

def plot_h2o_atomic_orbitals(centers):
    H1s = read_orbital_points("cpp/Hs000_0.000_-1.426_-0.888.txt")
    H2s = read_orbital_points("cpp/Hs000_0.000_1.426_-0.888.txt")
    O1px = read_orbital_points("cpp/Op100_0.000_0.000_0.222.txt")
    O1py = read_orbital_points("cpp/Op010_0.000_0.000_0.222.txt")
    O1pz = read_orbital_points("cpp/Op001_0.000_0.000_0.222.txt")
    O2s = read_orbital_points("cpp/Os000_0.000_0.000_0.222.txt")

    # Note the redundant p orbitals because of the two different signed isosurfaces
    # This is a hack of course: there is a cleaner factoring that I haven't come up with yet
    test_points = [H1s, H2s, O1px, O1py, O1pz, O2s]

    s_level = np.max(test_points[0]) * 0.5
    p_level = np.max(test_points[2]) * 0.4
    levels = [s_level, s_level, p_level, p_level, p_level, s_level]

    colors = ["blue", "blue", "red", "green", "yellow", "blue"]

    current_plot_development.plot_filtered_surfaces(test_points, levels, colors, centers, (-10.0, 10.0), True)


def plot_molecular_orbital(path_to_data, centers, default_color = "red"):
    print(path_to_data)
    MO_points = read_orbital_points(path_to_data)
    print(np.max(MO_points))
    levels = [np.max(MO_points) * 0.5]
    colors = [default_color]
    try:
        current_plot_development.plot_filtered_surfaces([MO_points], levels, colors, centers, (-10.0, 10.0), True)
    except:
        levels = [np.max(MO_points) * 200]
        current_plot_development.plot_filtered_surfaces([MO_points], levels, colors, centers, (-10.0, 10.0), True)


def plot_alpha_orbitals(directory, centers):
    # Search for files that match the given pattern in the directory
    pattern = "_alpha_MO_"
    files = []
    for filename in os.listdir(directory):
        if pattern in filename and filename.endswith(".txt"):
            files.append(filename)
    
    # Sort files based on the numerical part in the filename
    files_sorted = sorted(files)

    # Plot each file in sorted order
    for filename in files_sorted:
        full_path = os.path.join(directory, filename)
        plot_molecular_orbital(full_path, centers)

def plot_beta_orbitals(directory, centers):
    # Search for files that match the given pattern in the directory
    pattern = "_beta_MO_"
    files = []
    for filename in os.listdir(directory):
        if pattern in filename and filename.endswith(".txt"):
            files.append(filename)
    
    # Sort files based on the numerical part in the filename
    files_sorted = sorted(files)

    # Plot each file in sorted order
    for filename in files_sorted:
        full_path = os.path.join(directory, filename)
        plot_molecular_orbital(full_path, centers)


def plot_electron_density(directory, centers):
    pattern = "electron_densities"
    files = []
    for filename in os.listdir(directory):
        if pattern in filename and filename.endswith(".txt"):
            files.append(filename)

    print(len(files))
    assert(len(files) == 1), "Only expected one density file for this function"
    full_path = os.path.join(directory, files[0])

    density_array = read_orbital_points(full_path)
    levels = [np.max(density_array) * 0.4]
    colors = ["red"]
    current_plot_development.plot_filtered_surfaces([density_array], levels, colors, centers, (-10.0, 10.0), True)


if __name__ == "__main__":
    ho_centers = [(0.0, 0.0, 0.0), (1.830, 0.0, 0.0)] # Manually added for the black dots in the plots
    # plot_ho_atomic_orbitals(ho_centers)
    # plot_alpha_orbitals("cpp", ho_centers)
    # plot_electron_density("cpp", ho_centers)

    h2o_centers = [(0.0, 0.0, 0.0), (0.0, -1.426, -0.888), (0.0, 1.426, -0.888)]
    # plot_h2o_atomic_orbitals(h2o_centers)
    plot_alpha_orbitals("cpp", h2o_centers)
    # plot_electron_density("cpp", h2o_centers)