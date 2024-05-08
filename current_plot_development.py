import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import marching_cubes

POSITIVE_ISOSURFACE_OPACITY = 0.7
NEGATIVE_ISOSURFACE_OPACITY = 0.3

def add_isosurface(ax, data, level, color, alpha, grid_range=(-1, 1)):
    n_points = data.shape[0]
    linspace = np.linspace(grid_range[0], grid_range[1], n_points)

    verts, faces, normals, values = marching_cubes(data, level)
    scale_factor = (linspace[-1] - linspace[0]) / (n_points - 1)
    verts = verts * scale_factor + linspace[0]

    ax.plot_trisurf(verts[:, 1], verts[:, 2], verts[:, 0], triangles=faces, color=color, alpha=alpha)

    return np.max(np.abs(verts), axis=0)  # Return max absolute values of the vertices for range calculation

def plot_filtered_surfaces(datasets, levels, colors, centers, title, grid_range=(-1, 1), plot_negative_isos = False):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    max_range = np.array([0, 0, 0])

    for data, level, color in zip(datasets, levels, colors):
        max_verts = add_isosurface(ax, data, level, color, POSITIVE_ISOSURFACE_OPACITY, grid_range)
        max_range = np.maximum(max_range, max_verts)

        if plot_negative_isos:
            try:
                add_isosurface(ax, data, -level, color, NEGATIVE_ISOSURFACE_OPACITY, grid_range)  # Plot negative isosurface
            except Exception as e:
                pass 
                # print(f"Unable to plot negative isosurface for level {level}: {e}")
                # The above print was confusing, it implied a problem but  not all orbitals have negative isosurfaces


    # Plotting the center points after isosurface plotting to avoid depending on their count
    for center in centers:
        ax.scatter(*center, color='black', s=100)

    max_limit = max(max_range) * 1.1  # Adding 10% buffer to the maximum range found
    ax.set_xlim(-max_limit, max_limit)
    ax.set_ylim(-max_limit, max_limit)
    ax.set_zlim(-max_limit, max_limit)
    ax.set_aspect('equal')  # Ensure equal aspect ratio
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.title(title)
    plt.show()

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
    # Specify colors for each dataset
    colors = ['blue', 'green']
