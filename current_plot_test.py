import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import marching_cubes

def plot_filtered_surfaces(datasets, level, colors, centers, grid_range=(-1, 1)):
    # datasets is a list of 3D numpy arrays, each representing a different Gaussian field
    # colors is a list of colors corresponding to each dataset
    # centers is a list of tuples, each representing the center of the corresponding Gaussian dataset

    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for data, color, center in zip(datasets, colors, centers):
        n_points = data.shape[0]  # Handle varying sizes within the datasets
        linspace = np.linspace(grid_range[0], grid_range[1], n_points)

        verts, faces, normals, values = marching_cubes(data, level)

        scale_factor = (linspace[-1] - linspace[0]) / (n_points - 1)
        verts = verts * scale_factor + linspace[0]  # First scale and shift to the grid range
        offset = center - np.array([linspace[int(n_points/2)] for _ in range(3)]) * scale_factor
        verts += offset


        # Plotting each Gaussian's isosurface
        ax.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2], triangles=faces, color=color, alpha=0.6)

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.title('3D Isosurfaces of Multiple Gaussian Fields')

    plt.show()

def generate_3d_gaussian_data(center, sigma, grid_size=50):
    # Generate Gaussian data centered at 'center' with standard deviation 'sigma'
    x = np.linspace(-1, 1, grid_size)
    y = np.linspace(-1, 1, grid_size)
    z = np.linspace(-1, 1, grid_size)
    x, y, z = np.meshgrid(x, y, z)

    gaussian = np.exp(-(((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2) / (2. * sigma**2)))
    return gaussian

if __name__ == "__main__":
    # Define the centers and sigmas for Gaussian distributions
    centers = [(0, 0, 0), (0.5, 0.5, 0.5)]
    sigmas = [0.5, 0.3]

    # Generate the Gaussian datasets
    datasets = [generate_3d_gaussian_data(center, sigma, grid_size=50 + i * 10) for i, (center, sigma) in enumerate(zip(centers, sigmas))]

    # Define a common level to plot
    level = np.max(datasets[0]) * 0.5

    # Specify colors for each dataset
    colors = ['blue', 'green']

    # Plot the isosurfaces with specified colors
    plot_filtered_surfaces(datasets, level, colors, centers)