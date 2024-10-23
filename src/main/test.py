import os
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.Execution.BasicRunner import BasicRunner
from PyFoam.Execution.UtilityRunner import UtilityRunner
from func import dir, generate_mesh, build_case, save_inputs, surface_area, collect_results, plots_cost_pressure, pareto

import vtk
import numpy
import csv

def vtk_results_test(sn, vtk_dir, csv_file, area):
  """
  Store all vtk data (in vtk_dir) to a CSV file.
  """
  vtk_files = [f for f in os.listdir(vtk_dir) if f.endswith('.vtk')]
  if not vtk_files:
    print(f"No VTK files found in the directory: {vtk_dir}")
    return
  # Sort vtk_files based on the numeric part of the filename
  vtk_files.sort(key=lambda f: int(f.split('_')[-1].replace('.vtk', '')))
  with open(csv_file, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    # Write header only once
    csvwriter.writerow(['sn','time', 'x', 'y', 'z', 'p', 'Ux', 'Uy', 'Uz', 'area'])
  for vtk_file in vtk_files:
    vtk_path = os.path.join(vtk_dir, vtk_file)
    print(f"Reading VTK file '{vtk_path}'")
    # Create vtk reader object
    vtk_lector = vtk.vtkUnstructuredGridReader()
    vtk_lector.SetFileName(vtk_path)
    vtk_lector.Update()
    # Get VTK data
    vtk_data = vtk_lector.GetOutput()
    if not vtk_data:
      print("Failed to read VTK file.")
      return
    if vtk_data.IsA('vtkUnstructuredGrid'):
      # Extract points
      points = vtk_data.GetPoints()
      if not points:
        print("No points found in the VTK data.")
        return
      num_points = points.GetNumberOfPoints()
      points_array = numpy.zeros((num_points, 3))
      for i in range(num_points):
        points_array[i] = points.GetPoint(i)
      # Extract point data for pressure and velocity
      point_data = vtk_data.GetPointData()
      pressure_array = point_data.GetArray('p')
      velocity_array = point_data.GetArray('U')
      if not pressure_array:
        print("No pressure array found in the VTK data.")
        return
      if not velocity_array:
        print("No velocity array found in the VTK data.")
        return
      pressure_values = [pressure_array.GetValue(i) for i in range(num_points)]
      velocity_values = [velocity_array.GetTuple3(i) for i in range(num_points)]
      # Extract time from the file name
      time_step = vtk_file.split('_')[-1].replace('.vtk', '')
      # Append to the CSV file
      with open(csv_file, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write data rows
        for i in range(num_points):
          row = [
            sn,
            time_step,
            points_array[i][0], points_array[i][1], points_array[i][2],
            pressure_values[i],
            velocity_values[i][0], velocity_values[i][1], velocity_values[i][2],
            area
          ]
          csvwriter.writerow(row)
    else:
      print(f"Unsupported VTK data type in file: {vtk_file}")
  print(f"\nResults for simulation {sn} writen to '{csv_file}'")

def testrun(logfile=0):
  """
  Complete single simulation run for testing purposes.
  """
  pwd = os.path.abspath(os.getcwd()) # Get current working directory
  print("\n\033[36m[CREATING DIRECTORY]\033[0m\n")
  path = dir(pwd) # Set simulation directory
  case = SolutionDirectory(path) # Create foam case file
  os.chdir(path) # Change directory to simulation directory
    
  # Generate mesh
  print("\n\033[33m[MESHING]\033[0m\n")
  mesh_file = generate_mesh()

  # Generate base case files
  print("\n\033[92m[CONFIGURATION]\033[0m\n")
  build_case()
  save_inputs(f"{path}_params.json") # Save generated parameters to file

  # Convert mesh to foam
  print("\n\033[35m[CONVERTING MESH TO FOAM]\033[0m\n")
  msh_runner = BasicRunner(argv=["gmshToFoam", "-case", ".", mesh_file], silent=False)
  msh_runner.start()
  if msh_runner.runOK():
    print("\n\033[38;5;46m[✓]\033[0m Successful mesh conversion to foam.\n")
  else:
    os.chdir("..") # Exit simulation directory
    raise ValueError("\n\033[31m[X]\033[0m Unsuccessful mesh conversion to foam.\n") 

  # Run the simulation using icoFoam
  print("\n\033[31m[SIMULATION]\033[0m\n")
  if logfile:
    runner = BasicRunner(argv=["icoFoam", "-case", ".", ">", "log.icoFoam"], silent=False)
  else:
    runner = BasicRunner(argv=["icoFoam", "-case", "."], silent=False)
  runner.start()
  if runner.runOK():
    print("\n\033[38;5;46m[✓]\033[0m Simulation ran successfuly.\n")
  else:
    os.chdir("..") # Exit simulation directory
    raise ValueError("\n\033[31m[X]\033[0m Simulation failed. Check log files for more details.\n") 
  
  # Convert to VTK results
  print("\n\033[32m[CONVERTING FOAM TO VTK]\033[0m\n")
  vtk_runner = UtilityRunner(argv=["foamToVTK", "-case", "."])
  vtk_runner.start()
  if vtk_runner.runOK():
    print("\n\033[38;5;46m[✓]\033[0m Successful foam conversion to VTK.\n")
  else:
    os.chdir("..") # Exit simulation directory
    raise ValueError("\n\033[31m[X]\033[0m Unsuccessful foam conversion to VTK.\n") 
 
  # Process and store VTK results
  print("\n\033[34m[PROCESSING VTK]\033[0m\n")
  sim_data = "results1_sim_vtk.csv"
  collected_data = "results2_cost_avgp.csv"
  optimal_data = "results3_pareto.csv"
  area = surface_area()
  vtk_results_test(0, "VTK", sim_data, area)
  # Build CSV with area and average pressure data
  collect_results(sim_data,collected_data)
  plots_cost_pressure(collected_data)
  pareto(collected_data, optimal_data)

  os.chdir("..") # Exit simulation directory

testrun()
