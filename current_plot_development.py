import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skimage.measure import marching_cubes


def plot_filtered_surfaces(datasets, levels, colors, centers, grid_range=(-1, 1)):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    max_range = np.array([0, 0, 0])

    for data, level, color, center in zip(datasets, levels, colors, centers):
        n_points = data.shape[0]
        linspace = np.linspace(grid_range[0], grid_range[1], n_points)

        verts, faces, normals, values = marching_cubes(data, level)

        scale_factor = (linspace[-1] - linspace[0]) / (n_points - 1)
        verts = verts * scale_factor + linspace[0]

        # original_center = np.mean(verts, axis=0)
        # offset = np.array(center) - original_center
        # verts += offset

        for i in range(3):
            max_range[i] = max(max_range[i], np.abs(verts[:, i]).max())

        ax.plot_trisurf(verts[:, 1], verts[:, 0], verts[:, 2], triangles=faces, color=color, alpha=0.6)
        ax.scatter(*center, color='black', s=100)

    max_limit = max(max_range)
    ax.set_xlim(-max_limit, max_limit)
    ax.set_ylim(-max_limit, max_limit)
    ax.set_zlim(-max_limit, max_limit)
    ax.set_aspect('equal')
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

def plot_point_cloud(datasets, grid_range=(-1, 1), value_threshold=0.1):
    # This function plots each point in the dataset in a scatter plot where the size of the point
    # is proportional to the value at that point in the dataset.
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    for data in datasets:
        n_points = data.shape[0]
        linspace = np.linspace(grid_range[0], grid_range[1], n_points)
        
        # We create a mesh grid and flatten it to use in scatter
        x = np.linspace(grid_range[0], grid_range[1], n_points)
        y = np.linspace(grid_range[0], grid_range[1], n_points)
        z = np.linspace(grid_range[0], grid_range[1], n_points)
        x, y, z = np.meshgrid(x, y, z)
        
        # Flatten the arrays so they can be used in scatter plot
        x = x.flatten()
        y = y.flatten()
        z = z.flatten()
        values = data.flatten()

        # Filter points where values are below the threshold
        mask = values > value_threshold
        x = x[mask]
        y = y[mask]
        z = z[mask]
        values = values[mask]
        
        # Scale point sizes by some factor
        scaled_sizes = values / values.max() * 100  # Adjust 100 to scale sizes up or down
        
        ax.scatter(x, y, z, c='red', s=scaled_sizes, alpha=0.6)
        
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.title('Point Cloud Representation of Data')
    plt.show()


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