def plot_meshing_time_3D_surface(path, results=0):
  import json
  import numpy as np
  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import Axes3D
  if results == 0:
    with open(path, 'r') as json_file:
      data = json.load(json_file)
    min_mesh_sizes = [item['small_inlet_minimum_mesh_size'] for item in data]
    max_mesh_sizes = [item['small_inlet_maximum_mesh_size'] for item in data]
    times = [item['time'] for item in data]
  else:
    min_mesh_sizes = [item['small_inlet_minimum_mesh_size'] for item in results]
    max_mesh_sizes = [item['small_inlet_maximum_mesh_size'] for item in results]
    times = [item['time'] for item in results]
  # Create a sorted index based on time and then max_mesh_sizes
  sorted_indices = sorted(range(len(times)), key=lambda i: (times[i], max_mesh_sizes[i]))
  # Sort data based on the sorted indices
  times = [times[i] for i in sorted_indices]
  min_mesh_sizes = [min_mesh_sizes[i] for i in sorted_indices]
  max_mesh_sizes = [max_mesh_sizes[i] for i in sorted_indices]
  # Create 3D plot
  fig = plt.figure(figsize=(14, 10))
  ax = fig.add_subplot(111, projection='3d')
  # Plot the surface with a higher contrast colormap and antialiasing
  trisurf = ax.plot_trisurf(max_mesh_sizes, min_mesh_sizes, times, cmap='hsv', linewidth=0.5, edgecolor='k', antialiased=True)
  # Customize the viewing angle to better appreciate the surface
  #ax.view_init(elev=30, azim=120)  # Adjust elevation and azimuth angle
  # Add color bar to emphasize contrast in values
  fig.colorbar(trisurf, ax=ax, shrink=0.5, aspect=10)
  # Set labels and title
  ax.set_xlabel('Overall Mesh Size')
  ax.set_ylabel('Small Inlet Mesh Size')
  ax.set_zlabel('Meshing Time (s)')
  ax.set_title('3D Plot of Meshing Time vs Mesh Sizes')
  # Save and display plot
  #plt.savefig("meshing_time_plot_3d_surface2.pdf")
  plt.show()

import sys
if len(sys.argv) > 1:
  print("The first argument must be a path to a JSON file with mesh size and timming data. All other arguments provided will be ignored.")
  plot_meshing_time_3D_surface(sys.argv[1])
else:
  print("No argument was provided.")

