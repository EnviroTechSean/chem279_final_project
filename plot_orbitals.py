import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def generate_h2_points():
    return [(0, 0, 0)] * 41

"""
If dragging the scatter is no good, azimuth and elevation view inits
With like 3 different views, should be good
"""

# This technique shimmers! Not great
def print_random_scatter():
    # Example data
    x = np.random.rand(100) * 10  # X coordinates
    y = np.random.rand(100) * 10  # Y coordinates
    z = np.random.rand(100) * 10  # Z coordinates
    values = np.random.rand(100)  # Decimal values [0, 1] for opacity

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # The lighter opacities were shimmering if overlapping, seeing if that works ok
    alphas = 0.5 + 0.5 * values 
    scatter = ax.scatter(x, y, z, color='darkblue', alpha=alphas)

    # Show plot
    plt.show()


if __name__ == "__main__":
    print_random_scatter()