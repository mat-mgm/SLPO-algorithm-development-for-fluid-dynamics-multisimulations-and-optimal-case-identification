def test_monte_carlo_validity(path):
  import math
  import json
  import statistics
  import matplotlib.pyplot as plt
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
  #sorted_times = sorted(times)  # Sort the meshing times in ascending order
  #cumulative_means = [statistics.mean(sorted_times[:i+1]) for i in range(len(sorted_times))]
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

import sys
if len(sys.argv) > 1:
  print("The first argument must be a path to a JSON file with mesh size and timming data. All other arguments provided will be ignored.")
  test_monte_carlo_validity(sys.argv[1])
else:
  print("No argument was provided.")

