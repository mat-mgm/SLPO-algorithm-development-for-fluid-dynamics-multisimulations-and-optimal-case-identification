import math
import random
import os
import time
import gmsh
import json
import statistics
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def randfp(low, high, precision):
  """
  Generate a random floating point number of arbitrary precision (less or equal precision to the float data type).
  """
  factor = 10 ** precision
  return int(random.uniform(low, high) * factor) / factor

precision = 8
inputs = {}
# Geometric constrains
inputs["boundary_conditions"] = ["big-inlet", "big-outlet", "small-inlet", "walls", "fluid"]
inputs.update({"ox": 0, "oy": 0, "oz": 0})      # origin point
inputs["a0"] = round(2.0 * math.pi, precision)  # circle arc
inputs["a1"] = round(math.pi / 2.0, precision)  # torus turn angle
inputs["l1"] = 0.125                            # outlet pipe length
inputs["l2"] = inputs["l1"]                     # inlet pipe length
inputs["l3"] = inputs["l1"]                     # small inlet pipe length
inputs["r1"] = randfp(0.02, 0.045, precision)   # outlet pipe radius 20-45mm
inputs["r2"] = randfp(0.02, 0.045, precision)   # turn pipe radius (torus)
inputs["r3"] = inputs["l1"] / 2.0               # turn radius (torus)
inputs["r4"] = inputs["r2"]                     # inlet radius
inputs["r5"] = randfp(0.0025, 0.005, precision) # small inlet radius 2.5-5mm inside geometry
inputs["r6"] = inputs["r5"]                     # small inlet radius 2.5-5mm outside geometry
# Mesh settings
inputs["points_per_curve"] = 100
inputs["small_inlet_minimum_mesh_size"] = 0.0015 # Mesh size of small inlet
inputs["small_inlet_maximum_mesh_size"] = 0.005 # Mesh size everywhere else
inputs["small_inlet_minimum_distance"] = 0.0
inputs["small_inlet_maximum_distance"] = 0.005

def monte_carlo():
  """
  Refresh the randomized inputs (Montecarlo method).
  """
  global inputs
  inputs["r1"] = randfp(0.02, 0.045, precision)   # outlet pipe radius 20-45mm
  inputs["r2"] = randfp(0.02, 0.045, precision)   # turn pipe radius (torus)
  inputs["r4"] = inputs["r2"]                     # inlet radius
  inputs["r5"] = randfp(0.0025, 0.005, precision) # small inlet radius 2.5-5mm inside geometry
  inputs["r6"] = inputs["r5"]                     # small inlet radius 2.5-5mm outside geometry
  inputs["points_per_curve"] = 100
  inputs["small_inlet_minimum_mesh_size"] = random.uniform(0.0025, 0.00025) # Mesh size of small inlet (0.0015)
  inputs["small_inlet_maximum_mesh_size"] = random.uniform(0.0075, 0.00075) # Mesh size everywhere else (0.005)
  inputs["big_ratio"] = inputs["r2"] / inputs["small_inlet_maximum_mesh_size"]
  inputs["small_ratio"] = inputs["r5"] / inputs["small_inlet_minimum_mesh_size"]
  inputs["small_inlet_minimum_distance"] = 0.0
  inputs["small_inlet_maximum_distance"] = 0.005

def generate_mesh():
  """
  Generate a mesh from input parameters and save to file.
  """
  gmsh.initialize()
  gmsh.model.add("fluid") # Create model
  # Create geometry
  if inputs["r1"] == inputs["r2"]:
    volume1 = gmsh.model.occ.addCylinder(inputs["ox"], inputs["oy"], inputs["oz"], 0, 0, inputs["l1"], inputs["r1"])
  else:
    volume1 = gmsh.model.occ.addCone(inputs["ox"], inputs["oy"], inputs["oz"], 0, 0, inputs["l1"], inputs["r1"], inputs["r2"])
  volume2 = gmsh.model.occ.addTorus(inputs["ox"] - inputs["r3"], inputs["oy"], inputs["oz"] + inputs["l1"], inputs["r3"], inputs["r2"], -1, inputs["a1"], zAxis=[0, 1, 0])
  if inputs["r2"] == inputs["r4"]:
    volume3 = gmsh.model.occ.addCylinder(inputs["ox"] + inputs["r3"] * (math.cos(inputs["a1"]) - 1), inputs["oy"], inputs["oz"] + inputs["l1"] + inputs["r3"] * math.sin(inputs["a1"]), -inputs["l2"] * math.sin(inputs["a1"]), 0, inputs["l2"] * math.cos(inputs["a1"]), inputs["r4"])
  else:
    volume3 = gmsh.model.occ.addCone(inputs["ox"] + inputs["r3"] * (math.cos(inputs["a1"]) - 1), inputs["oy"], inputs["oz"] + inputs["l1"] + inputs["r3"] * math.sin(inputs["a1"]), -inputs["l2"] * math.sin(inputs["a1"]), 0, inputs["l2"] * math.cos(inputs["a1"]), inputs["r2"], inputs["r4"])
  if inputs["r5"] == inputs["r6"]:
    volume4 = gmsh.model.occ.addCylinder(inputs["ox"] + inputs["r3"] * (math.cos(inputs["a1"] / 2) - 1), inputs["oy"], inputs["oz"] + inputs["l1"] + inputs["r3"] * math.sin(inputs["a1"] / 2), inputs["l3"] * math.cos(inputs["a1"]), 0, inputs["l3"] * math.sin(inputs["a1"]), inputs["r5"])
  else:
    volume4 = gmsh.model.occ.addCone(inputs["ox"] + inputs["r3"] * (math.cos(inputs["a1"] / 2) - 1), inputs["oy"], inputs["oz"] + inputs["l1"] + inputs["r3"] * math.sin(inputs["a1"] / 2), inputs["l3"] * math.cos(inputs["a1"]), 0, inputs["l3"] * math.sin(inputs["a1"]), inputs["r5"], inputs["r6"])
  gmsh.model.occ.synchronize() # Sync changes
  volumes = [volume1, volume2, volume3, volume4]
  volume_tuples = [(3, v) for v in volumes]
  fused_volume, _ = gmsh.model.occ.fuse(volume_tuples, volume_tuples)
  gmsh.model.occ.synchronize() # Sync changes
  ## Physical boundaries
  gmsh.model.addPhysicalGroup(2, [7], 14)
  gmsh.model.setPhysicalName(2, 14, inputs["boundary_conditions"][0])
  gmsh.model.addPhysicalGroup(2, [3], 15)
  gmsh.model.setPhysicalName(2, 15, inputs["boundary_conditions"][1])
  gmsh.model.addPhysicalGroup(2, [6], 16)
  gmsh.model.setPhysicalName(2, 16, inputs["boundary_conditions"][2])
  gmsh.model.addPhysicalGroup(2, [2, 5, 1, 4], 17)
  gmsh.model.setPhysicalName(2, 17, inputs["boundary_conditions"][3])
  gmsh.model.addPhysicalGroup(3, [fused_volume[0][1]], 18)
  gmsh.model.setPhysicalName(3, 18, inputs["boundary_conditions"][4]) # fluid volume
  gmsh.model.occ.synchronize() # Sync changes
  # Define mesh size fields for the small inlet
  gmsh.model.mesh.field.add("Distance", 1)
  gmsh.model.mesh.field.setNumbers(1, "FacesList", [4, 6])
  gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", inputs["points_per_curve"])
  gmsh.model.mesh.field.add("Threshold", 2)
  gmsh.model.mesh.field.setNumber(2, "InField", 1)
  gmsh.model.mesh.field.setNumber(2, "SizeMin", inputs["small_inlet_minimum_mesh_size"])
  gmsh.model.mesh.field.setNumber(2, "SizeMax", inputs["small_inlet_maximum_mesh_size"])
  gmsh.model.mesh.field.setNumber(2, "DistMin", inputs["small_inlet_minimum_distance"])
  gmsh.model.mesh.field.setNumber(2, "DistMax", inputs["small_inlet_maximum_distance"])
  gmsh.model.mesh.field.setAsBackgroundMesh(2)
  # Set options
  gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
  gmsh.option.setNumber("General.Terminal", 1) # Enable terminal output
  gmsh.option.setNumber("General.NumThreads", os.cpu_count()) # Use all threads for parallel processing
  if inputs["mesh_type"]: # quad
    gmsh.option.setNumber("Mesh.Algorithm", 6) # Use Front-Delaunay algorithm
    gmsh.option.setNumber("Mesh.RecombineAll", 1) # Enable recombination of triangular elements into quadrilaterals
    surfaces = gmsh.model.getEntities(dim=2)
    for surface in surfaces:
      gmsh.model.mesh.setRecombine(2, surface[1])
      gmsh.model.mesh.setSmoothing(2, surface[1], 10) # Set smoothing to improve element quality
  gmsh.model.mesh.generate(3) # Generate 3D mesh
  gmsh.finalize()

def profile_machine(path):
  import json
  import platform
  import psutil
  import subprocess
  info = {}
  # CPU information
  cpu_info = {}
  cpu_info['brand'] = platform.processor()
  cpu_info['count'] = psutil.cpu_count(logical=False)
  cpu_info['threads'] = psutil.cpu_count(logical=True)
  cpu_info['frequency'] = psutil.cpu_freq().max
  info['cpu'] = cpu_info
  # Memory information
  mem_info = psutil.virtual_memory()
  info['memory'] = {
    'total_gb': mem_info.total / (1024 ** 3)
  }
  # Motherboard information
  try:
    motherboard_info = subprocess.check_output(['dmidecode', '-t', 'baseboard']).decode()
    motherboard_brand = ''
    motherboard_model = ''
    for line in motherboard_info.split('\n'):
      if 'Manufacturer:' in line:
        motherboard_brand = line.split(':')[1].strip()
      elif 'Product Name:' in line:
        motherboard_model = line.split(':')[1].strip()
    info['motherboard'] = {
      'brand': motherboard_brand,
      'model': motherboard_model
    }
  except Exception as e:
    info['motherboard'] = 'dmidecode command not available or not permitted'
  # System information
  system_info = platform.uname()
  info['system'] = {
    'system': system_info.system,
    'hostname': system_info.node,
    'release': system_info.release,
    'version': system_info.version,
    'architecture': system_info.machine,
    'processor': system_info.processor
  }
  # Save to file
  with open(path, 'w') as f:
    json.dump(info, f, indent=2)
  print("Machine profile generated.")

def plot_meshing_time_3D_surface(path, results=0):
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
  # Create a meshgrid for the surface plot
  min_mesh_sizes_unique = np.unique(min_mesh_sizes)
  max_mesh_sizes_unique = np.unique(max_mesh_sizes)
  X, Y = np.meshgrid(min_mesh_sizes_unique, max_mesh_sizes_unique)
  # Interpolate Z values on the grid
  Z = np.array([[np.interp(x, min_mesh_sizes, times) for x in min_mesh_sizes_unique] for y in max_mesh_sizes_unique])
  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111, projection='3d')
  # Plot the surface
  ax.plot_surface(X, Y, Z, cmap='viridis')
  # Reverse the direction of the small inlet mesh size axis
  ax.set_xlim(max(min_mesh_sizes), min(min_mesh_sizes))
  # Set axis labels
  ax.set_xlabel('Small Inlet Mesh Size')
  ax.set_ylabel('Overall Mesh Size')
  ax.set_zlabel('Time (s)')
  ax.set_title('3D Surface Plot of Mesh Sizes vs. Time')
  plt.savefig("meshing_time_plot_3d_surface.pdf")

def test_monte_carlo_validity(path):
  try:
    with open(path, "r") as json_file:
      data = json.load(json_file)
  except Exception as e:
    print(f"Error loading {path}: {e}")
    return
  # Extract relevant data for analysis
  times = [entry["time"] for entry in data if "time" in entry]
  r1_values = [entry["r1"] for entry in data if "r1" in entry]
  r2_values = [entry["r2"] for entry in data if "r2" in entry]
  r4_values = [entry["r4"] for entry in data if "r4" in entry]
  r5_values = [entry["r5"] for entry in data if "r5" in entry]
  r6_values = [entry["r6"] for entry in data if "r6" in entry]
  small_inlet_minimum_mesh_size = [entry["small_inlet_minimum_mesh_size"] for entry in data if "small_inlet_minimum_mesh_size" in entry]
  small_inlet_maximum_mesh_size = [entry["small_inlet_maximum_mesh_size"] for entry in data if "small_inlet_maximum_mesh_size" in entry]
  if not times:
    print("No valid 'time' data found in the JSON file.")
    return
  # Perform statistical analysis
  mean_time = statistics.mean(times)
  stdev_time = statistics.stdev(times) if len(times) > 1 else 0
  min_time = min(times)
  max_time = max(times)
  # Print the results
  print(f"Number of samples: {len(times)}")
  print(f"Mean meshing time: {mean_time:.2f} seconds")
  print(f"Standard deviation of meshing time: {stdev_time:.2f} seconds")
  print(f"Minimum meshing time: {min_time:.2f} seconds")
  print(f"Maximum meshing time: {max_time:.2f} seconds")
  # Example validity check: Mean and standard deviation should be within expected ranges
  expected_mean_range = (1.0, 10.0)  # Adjust as necessary
  expected_stdev_range = (0.0, 5.0)  # Adjust as necessary
  if not (expected_mean_range[0] <= mean_time <= expected_mean_range[1]):
    print(f"Mean meshing time is out of the expected range {expected_mean_range}.")
  else:
    print("Mean meshing time is within the expected range.")
  if not (expected_stdev_range[0] <= stdev_time <= expected_stdev_range[1]):
    print(f"Standard deviation of meshing time is out of the expected range {expected_stdev_range}.")
  else:
    print("Standard deviation of meshing time is within the expected range.")
  # Plotting to visualize distribution and convergence
  plt.figure(figsize=(14, 6))
  # Distribution of meshing times
  plt.subplot(1, 2, 1)
  plt.hist(times, bins=20, alpha=0.7, color='blue', edgecolor='black')
  plt.title('Distribution of Meshing Times')
  plt.xlabel('Meshing Time (seconds)')
  plt.ylabel('Frequency')
  # Convergence check: Plot mean meshing time over iterations
  cumulative_means = [statistics.mean(times[:i+1]) for i in range(len(times))]
  plt.subplot(1, 2, 2)
  plt.plot(cumulative_means, marker='o', linestyle='-')
  plt.title('Convergence of Mean Meshing Time')
  plt.xlabel('Number of Tests')
  plt.ylabel('Cumulative Mean Meshing Time (seconds)')
  plt.tight_layout()
  plt.savefig("statistical_plot.pdf")
  # Coverage check for input parameters
  def check_coverage(param_values, param_name, param_range):
    min_val, max_val = min(param_values), max(param_values)
    if min_val < param_range[0] or max_val > param_range[1]:
      print(f"Parameter {param_name} is out of the expected range {param_range}.")
    else:
      print(f"Parameter {param_name} covers the expected range {param_range}.")
  check_coverage(r1_values, "r1", (0.02, 0.045))
  check_coverage(r2_values, "r2", (0.02, 0.045))
  check_coverage(r4_values, "r4", (0.02, 0.045))
  check_coverage(r5_values, "r5", (0.0025, 0.005))
  check_coverage(r6_values, "r6", (0.0025, 0.005))
  check_coverage(small_inlet_minimum_mesh_size, "small_inlet_minimum_mesh_size", (0.00025, 0.0025))
  check_coverage(small_inlet_maximum_mesh_size, "small_inlet_maximum_mesh_size", (0.00075, 0.0075))

def dir(base_dir, prefix='s'):
    counter = 0
    while True:
        dir_name = f"{prefix}{counter}"
        dir_path = os.path.join(base_dir, dir_name)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(f"Created directory:  '{dir_name}'")
            return dir_name
        else:
            counter += 1

def read_json(file_path):
  """
  Function to read the existing JSON data
  """
  if os.path.exists(file_path):
    with open(file_path, 'r') as file:
      return json.load(file)
  return []
def write_json(file_path, data):
  """
  Function to write updated JSON data
  """
  with open(file_path, 'w') as file:
    json.dump(data, file, indent=2)
def append_to_json(file_path, new_data):
  """
  Function to append data to the JSON file
  """
  # Read existing data
  data = read_json(file_path)
  # Append new data
  data.append(new_data)
  # Write updated data
  write_json(file_path, data)

def plot_existing_data(path):
  plot_meshing_time_3D_surface(path)
  test_monte_carlo_validity(path)

def mesh_test():
  global inputs
  inputs["mesh_type"] = 0
  number_of_tests = 100
  dirname = dir(".", f"results_{inputs['mesh_type']}_n{number_of_tests}_")
  os.chdir(dirname) # Change to test directory
  print(f"Running {number_of_tests} meshing tests...")
  results = []
  path = "meshing_time.json"
  for i in range(number_of_tests):
    try:
      monte_carlo()
      t0 = time.time() # Set time 0
      generate_mesh() # Generate mesh
      inputs["time"] = time.time() - t0
      append_to_json(path, inputs.copy())
      print(f"{i+1} out of {number_of_tests}")
    except KeyboardInterrupt:
      print(f"Last mesh size was:\nMesh size: {inputs['small_inlet_maximum_mesh_size']}\nSmall inlet mesh size: {inputs['small_inlet_minimum_mesh_size']}")
  plot_existing_data(path)
  os.chdir("..")
  print("Test finished")

mesh_test()
#plot_existing_data("meshing_time_tri.json")
