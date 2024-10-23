def plot_meshing_time_2d(path, results=0):
  import json
  import matplotlib.pyplot as plt
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
  plt.figure(figsize=(10, 5))
  plt.subplot(1, 2, 1)
  plt.scatter(min_mesh_sizes, times, color='blue', label='Max Mesh Size')
  plt.xlabel('Maximum Mesh Size')
  plt.ylabel('Time (s)')
  plt.title('Small Inlet Mesh Size vs. Time')
  plt.legend()
  plt.subplot(1, 2, 2)
  plt.scatter(max_mesh_sizes, times, color='red', label='Max Mesh Size')
  plt.xlabel('Maximum Mesh Size')
  plt.ylabel('Time (s)')
  plt.title('Overall Mesh Size vs. Time')
  plt.legend()
  plt.tight_layout()
  plt.savefig("meshing_time_plot.pdf")
  #plt.show()

import sys
if len(sys.argv) > 1:
  print("The first argument must be a path to a JSON file with mesh size and timming data. All other arguments provided will be ignored.")
  plot_meshing_time_2d(sys.argv[1])
else:
  print("No argument was provided.")
