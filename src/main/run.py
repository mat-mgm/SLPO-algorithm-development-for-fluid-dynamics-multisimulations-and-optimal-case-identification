import os
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.Execution.BasicRunner import BasicRunner
from PyFoam.Execution.UtilityRunner import UtilityRunner
from func import dir, cpdir, generate_mesh, build_case, mod_case, monte_carlo, surface_area, save_inputs, csv_template, result_region, vtk_results, collect_results, pareto, plots_cost_pressure, compute_statistics

# Filenames for templates and results
sim_data = "results1_sim_vtk.csv"
collected_data = "results2_avgp_cost.csv"
optimal_data = "results3_pareto.csv"

def initrun():
  """
  Initial run without mesh generation or simulation to build the template case and the results file.
  """
  print("\n\033[34m[INITIAL RUN ➡ ]\033[0m\nTemplate case generation:")
  pwd = os.path.abspath(os.getcwd()) # Get current working directory
  print("\n\033[36m[CREATING DIRECTORY]\033[0m\n")
  path = dir(pwd) # Set simulation directory
  os.chdir(path) # Change directory to simulation directory

  # Generate base case files
  print("\n\033[92m[CONFIGURATION]\033[0m\n")
  build_case()
  os.chdir("..") # Exit simulation directory

  # Generate CSV files for storing results
  csv_template(sim_data, ['sn','time', 'x', 'y', 'z', 'p', 'Ux', 'Uy', 'Uz','area'])
  csv_template(collected_data, ['sn','avg_p','area'])
  
  return path

def stdrun(template, number):
  """
  Standard run which copies the case files from an initial run and runs a full simulation.
  """
  print(f"\n\033[34m[RUN {number} ➡ ]\033[0m\n")
  pwd = os.path.abspath(os.getcwd()) # Get current working directory
  # Copy template files
  print("\n\033[36m[COPYING TEMPLATE]\033[0m\n")
  path = cpdir(pwd, template) # Copy template
  case = SolutionDirectory(path) # Create foam case file
  os.chdir(path) # Change directory to simulation directory

  print("\n\033[32m[MONTE CARLO]\032[0m\n")
  monte_carlo() # Generate new random values for the simulation
  save_inputs(f"{path}_params.json") # Save generated parameters to file

  # Generate mesh
  print("\n\033[33m[MESHING]\033[0m\n")
  mesh_file = generate_mesh()

  # Generate base case files
  print("\n\033[92m[CONFIGURATION]\033[0m\n")
  mod_case(".") # Add generated parameters to template case

  # Convert mesh to foam
  print("\n\033[35m[CONVERTING MESH TO FOAM]\033[0m\n")
  msh_runner = BasicRunner(argv=["gmshToFoam", "-case", ".", mesh_file], silent=False)
  msh_runner.start()
  if msh_runner.runOK():
    print("\n\033[38;5;46m[✓]\033[0m Successful mesh conversion to foam.\n")
  else:
    print("\n\033[31m[X]\033[0m Unsuccessful mesh conversion to foam.\n") 
    os.chdir("..") # Exit simulation directory
    return

  # Run the simulation using icoFoam
  print("\n\033[31m[SIMULATION]\033[0m\n")
  runner = BasicRunner(argv=["icoFoam", "-case", "."], silent=False)
  runner.start()
  if runner.runOK():
    print("\n\033[38;5;46m[✓]\033[0m Simulation ran successfuly.\n")
  else:
    print("\n\033[31m[X]\033[0m Simulation failed. Check log files for more details.\n") 
    os.chdir("..") # Exit simulation directory
    return

  # Convert to VTK results
  print("\n\033[32m[CONVERTING FOAM TO VTK]\033[0m\n")
  vtk_runner = UtilityRunner(argv=["foamToVTK", "-case", "."])
  vtk_runner.start()
  if vtk_runner.runOK():
    print("\n\033[38;5;46m[✓]\033[0m Successful foam conversion to VTK.\n")
  else:
    print("\n\033[31m[X]\033[0m Unsuccessful foam conversion to VTK.\n") 
    os.chdir("..") # Exit simulation directory
    return
 
  # Process and store VTK results in a CSV file
  print("\n\033[34m[PROCESSING VTK]\033[0m\n")
  region = result_region() # Restric results to points near th small inlet entrance
  area = surface_area()
  # Extract all VTK data from the specified geometrical region
  vtk_results(number, "VTK", os.path.join("..", sim_data), region[0], region[1], region[2], area)

  os.chdir("..") # Exit simulation directory

def main(i):
  """
  Flow control function to make an initial run generating a template case and arbitrary number i of consequent standard runs.
  """
  template = initrun()
  i+=1
  for num in range(1,i):
    stdrun(template, num)
  collect_results(sim_data, collected_data)
  plots_cost_pressure(collected_data)
  compute_statistics(collected_data)
  pareto(collected_data, optimal_data)

if __name__ == '__main__':
  main(2) # Number of simulations to perform

# To regenerate CSV and plots when simulations are already done
#collect_results(sim_data, collected_data)
#plots_cost_pressure(collected_data)
#compute_statistics(collected_data)
#pareto(collected_data,optimal_data)
