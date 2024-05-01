import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata

import current_plot_test

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