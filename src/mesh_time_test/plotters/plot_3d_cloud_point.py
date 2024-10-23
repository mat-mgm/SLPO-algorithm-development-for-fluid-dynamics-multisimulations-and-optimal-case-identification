def plot_meshing_time_3D_point_cloud(path, results=0):
  import json
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
  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(min_mesh_sizes, max_mesh_sizes, times, color='blue', label='Mesh Sizes vs Time')
  ax.set_xlabel('Small Inlet Mesh Size')
  ax.set_ylabel('Overall Mesh Size')
  ax.set_zlabel('Time (s)')
  ax.set_title('3D Plot of Mesh Sizes vs. Time')
  ax.legend()
  plt.savefig("meshing_time_plot_3d_point_cloud.pdf")

import sys
if len(sys.argv) > 1:
  print("The first argument must be a path to a JSON file with mesh size and timming data. All other arguments provided will be ignored.")
  plot_meshing_time_3D_point_cloud(sys.argv[1])
else:
  print("No argument was provided.")

