import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

"""
If dragging the scatter is no good, azimuth and elevation view inits
With like 3 different views, should be good
"""

# This technique shimmers! Not great
def plot_random_scatter():
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


def plot_meshgrid():
    # Creating the grid points
    x = np.linspace(-1, 1, 10)
    y = np.linspace(-1, 1, 10)
    z = np.linspace(-1, 1, 10)
    x, y, z = np.meshgrid(x, y, z)

    print(x.shape)

    # We'll only sample every 5th point from each dimension for clarity in visualization
    sample_rate = 5
    x_sample = x[::sample_rate, ::sample_rate, ::sample_rate]
    y_sample = y[::sample_rate, ::sample_rate, ::sample_rate]
    z_sample = z[::sample_rate, ::sample_rate, ::sample_rate]

    # Plotting
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x.flatten(), y.flatten(), z.flatten(), c='blue', alpha=0.6)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('3D Grid Points Visualization')
    plt.show()


def plot_gaussian_interpolation_slice():
    # Create a 3D grid
    x = np.linspace(-1, 1, 30)
    y = np.linspace(-1, 1, 30)
    z = np.linspace(-1, 1, 30)
    x, y, z = np.meshgrid(x, y, z)
    
    # Define a Gaussian distribution centered at 0,0,0
    sigma = 0.5
    gaussian = np.exp(-((x**2 + y**2 + z**2) / (2. * sigma**2)))
    
    # Prepare data for interpolation
    points = np.array([x.flatten(), y.flatten(), z.flatten()]).T
    values = gaussian.flatten()
    
    # Create a finer grid for interpolation
    grid_x, grid_y, grid_z = np.mgrid[-1:1:10j, -1:1:10j, -1:1:10j]
    
    # Interpolate using griddata
    interpolated_values = griddata(points, values, (grid_x, grid_y, grid_z), method='linear')

    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot a surface
    # We will slice the volume to extract a surface to plot, you can adjust this slice
    x_slice = grid_x[:, :, 9]  # Adjust the slice index as needed
    print(grid_x.shape)
    y_slice = grid_y[:, :, 9]  # Adjust the slice index as needed
    z_slice = interpolated_values[:, :, 9]  # Surface based on interpolated values

    # Masking NaN values to avoid issues in plotting
    mask = ~np.isnan(z_slice)
    x_slice = x_slice[mask]
    y_slice = y_slice[mask]
    z_slice = z_slice[mask]
    print(x_slice.shape)

    # Plot surface
    surf = ax.plot_trisurf(x_slice, y_slice, z_slice, cmap='viridis', linewidth=0.1)
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    plt.title('3D Gaussian Surface Slice')
    plt.show()


if __name__ == "__main__":
    # plot_random_scatter()
    # plot_meshgrid()
    # plot_gaussian_interpolation_slice()
    pass