import os
import os.path
import re
import sys

import subprocess # subprocess is for calling the c++ calculations, see generate_molecule_points

import numpy as np

import current_plot_development

def call_mo_points_main(molecule_file : str):
    executable_path = "./cpp/mo_points_main"
    assert os.path.exists(executable_path), "Expected this script to run from top level directory, with 'make all' having been run from cpp directory beforehand."
    assert os.path.exists(molecule_file), f"Input molecule file {molecule_file} not found, run from top level directory"

    cmd = [executable_path, molecule_file, "./basis"]
    
    try:
        result = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        print("Error occurred while executing the command:", e)
        return e.output

def make_clean():
    assert os.path.exists("./Makefile"), "Expected Makefile with clean target"

    cmd = ["make", "clean"]  # Corrected command list
    try:
        result = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        print("Error occurred while executing the command:", e)
        return e.output

def read_orbital_points(filepath):
    data = np.loadtxt(filepath)
    data = data.reshape(50, 50, 50)
    return data


def extract_centers(directory):
    pattern = r"_(-?\d+\.\d{3})_(-?\d+\.\d{3})_(-?\d+\.\d{3})\.txt"

    centers = set()

    for filename in os.listdir(directory):
        match = re.search(pattern, filename)
        if match:
            centers.add(tuple(map(float, match.groups())))

    return list(centers)

def plot_atomic_orbitals(directory, centers):
    color_map = {"s": "blue", "p": "red", "d": "green", "f": "yellow"}
    level_map = {"s": 0.6, "p": 0.4, "d": 0.3, "f": 0.2}

    test_points = []
    colors = []
    levels = []

    print(len(centers))

    for center in centers:
        for filename in os.listdir(directory):
            if all(f"{c:.3f}" in filename for c in center):  # Checking if filename contains the center coordinates
                filepath = os.path.join(directory, filename)
                points = read_orbital_points(filepath)
                test_points.append(points)
                max_point = np.max(points)
                orbital_type = filepath.split('/')[-1][1]  # Extracting the orbital type from filename
                level = max_point * level_map.get(orbital_type, 0.5)  # Using a default level if type not found
                colors.append(color_map.get(orbital_type, "gray"))  # Default color if not specified
                levels.append(level)

    current_plot_development.plot_filtered_surfaces(test_points, levels, colors, centers, "Atomic Orbitals", (-10.0, 10.0), True)

def plot_molecular_orbitals(directory, centers):
    pattern = "alpha_MO"
    files = [f for f in os.listdir(directory) if pattern in f]
    for filename in files:
        filepath = os.path.join(directory, filename)
        points = read_orbital_points(filepath)
        max_point = np.max(points)
        level = max_point * 0.5
        color = "red"
        try:
            current_plot_development.plot_filtered_surfaces([points], [level], color, centers, "Molecular Orbital", (-10.0, 10.0), True)
        except: # Pretty trial-and-error here
            level = np.max(points) * 200
            current_plot_development.plot_filtered_surfaces([points], [level], color, centers, "Molecular Orbital", (-10.0, 10.0), True)

def plot_electron_density(directory, centers):
    pattern = "electron_densities"
    for filename in os.listdir(directory):
        if pattern in filename:
            filepath = os.path.join(directory, filename)
            points = read_orbital_points(filepath)
            max_point = np.max(points)
            level = max_point * 0.4
            color = "blue"
            current_plot_development.plot_filtered_surfaces([points], [level], [color], centers, "Electron Density", (-10.0, 10.0), True)
            break  # Assume only one electron density file

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_molecule_file>")
        sys.exit(1)

    molecule_file = sys.argv[1]
    directory = "./"  

    # Step 1: Clean previous files and generate new data
    make_clean()
    call_mo_points_main(molecule_file)

    # Step 2: Extract centers
    centers = extract_centers(directory)

    # Step 3: Plot atomic orbitals
    plot_atomic_orbitals(directory, centers)

    # Step 4: Plot Molecular Orbitals
    plot_molecular_orbitals(directory, centers)

    # Step 5: Plot Electron Densities
    plot_electron_density(directory, centers)