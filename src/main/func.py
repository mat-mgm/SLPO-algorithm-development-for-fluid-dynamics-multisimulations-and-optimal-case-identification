# FUNCTIONS FILE

# MODULES

import math
import random
import os
import json
import gmsh
import shutil
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import vtk
import numpy
import csv
import pandas as pd
import matplotlib.pyplot as plt
import statistics
from datetime import datetime

# MATHS

def randfp(low, high, precision):
  """
  Generate a random floating point number of arbitrary precision (less or equal precision to the float data type).
  """
  factor = 10 ** precision
  return int(random.uniform(low, high) * factor) / factor

def surface_area():
  """
  Calculate the surface area of the mesh based on the geometrical parameters from inputs dictionary.
  """
  global inputs
  area_components = []
  # Outlet
  if inputs["r1"] == inputs["r2"]: # lateral surface of cylinder A = 2πr(h+r)
      area_components.append(2 * math.pi * inputs["r1"] * (inputs["l1"] + inputs["r1"]))
  else: # lateral surface of conical frustum A = π(r1 + r2)√((r1 - r2)² + h²) 
      area_components.append(math.pi * (inputs["r1"] + inputs["r2"]) * math.sqrt((inputs["r1"] - inputs["r2"])**2 + inputs["l1"]**2))
  # Turn (torus)
  # surface of torus A = (2πr)(2πR)
  area_components.append((2 * math.pi * inputs["r2"]) * (inputs["a1"] * inputs["r3"]))
  # Inlet
  if inputs["r2"] == inputs["r4"]: # lateral surface of cylinder A = 2πr(h+r)
      area_components.append(2 * math.pi * inputs["r2"] * (inputs["l2"] + inputs["r2"]))
  else: # lateral surface of conical frustum A = π(r1 + r2)√((r1 - r2)² + h²)
      area_components.append(math.pi * (inputs["r2"] + inputs["r4"]) * math.sqrt((inputs["r2"] - inputs["r4"])**2 + inputs["l2"]**2))
  # Small inlet
  if inputs["r5"] == inputs["r6"]: # lateral surface of cylinder A = 2πr(h+r)
      area_components.append(2 * math.pi * inputs["r5"] * (inputs["l3"] + inputs["r5"]))
  else: # lateral surface of conical frustum A = π(r1 + r2)√((r1 - r2)² + h²)
      area_components.append(math.pi * (inputs["r5"] + inputs["r6"]) * math.sqrt((inputs["r5"] - inputs["r6"])**2 + inputs["l3"]**2))
  if not area_components:
    raise ValueError("No area components were calculated.")
  return sum(area_components)

def result_region():
  """
  Define the geometric range from which to extract simulation results.
  """
  z_offset = 0.001
  return [
    (inputs["ox"] + inputs["r3"] * (math.cos(inputs["a1"] / 2) - 1) + inputs["l3"] * math.cos(inputs["a1"]) - inputs["r6"], inputs["ox"] + inputs["r3"] * (math.cos(inputs["a1"] / 2) - 1) + inputs["l3"] * math.cos(inputs["a1"]) + inputs["r6"]),
    (inputs["oy"] - inputs["r6"], inputs["oy"] + inputs["r6"]),
    (inputs["oz"] + inputs["l1"] + inputs["r3"] * math.sin(inputs["a1"] / 2) + inputs["l3"] * math.sin(inputs["a1"]) - z_offset, inputs["oz"] + inputs["l1"] + inputs["r3"] * math.sin(inputs["a1"] / 2) + inputs["l3"] * math.sin(inputs["a1"]))
  ]

# INPUTS DICTIONARY

precision = 8
inputs = {}

# Simulation constrains
inputs["starttime"] = 0
inputs["endtime"] = 0.01
inputs["deltat"] =  0.0001
inputs["writeControl"] = "timeStep" # "runTime" or "timeStep"
inputs["writeinterval"] = 5 #inputs["deltat"] # integer 1 or bigger for timeInterval

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
inputs["area"] = surface_area() # Geometry external area

# Mesh settings
inputs["mesh_type"] = 0 # 0 for tri; 1 for quad
inputs["points_per_curve"] = 100
inputs["small_inlet_minimum_mesh_size"] = 0.001 # Mesh size of small inlet
if not inputs["mesh_type"]: # tri
  inputs["small_inlet_maximum_mesh_size"] = 0.005 # Mesh size everywhere else (equivalent to mesh_maximum_length)
  inputs["mesh_maximum_length"] = 0 # 0 means disabled
else: # quad
  inputs["small_inlet_maximum_mesh_size"] = 0.005
  inputs["mesh_maximum_length"] = 0.02
inputs["small_inlet_minimum_distance"] = 0.0
inputs["small_inlet_maximum_distance"] = 0.005

# Pysical parameters
p_atm = 101325 # Atmosferic pressure in pascals
def monte_carlo():
  """
  Refresh the randomized inputs (Monte Carlo method).
  """
  global inputs
  inputs["r1"] = randfp(0.02, 0.045, precision)   # outlet pipe radius 20-45mm
  inputs["r2"] = randfp(0.02, 0.045, precision)   # turn pipe radius (torus)
  inputs["r4"] = inputs["r2"]                     # inlet radius
  inputs["r5"] = randfp(0.0025, 0.005, precision) # small inlet radius 2.5-5mm inside geometry
  inputs["r6"] = inputs["r5"]                     # small inlet radius 2.5-5mm outside geometry
  inputs["area"] = surface_area() # Geometry external area
  inputs["temperature"] = randfp(268.15, 313.15, precision) # T [K] (-5ºC, 40ºC)
  inputs["big_outlet_presure"] = randfp(p_atm, p_atm-300, precision)
  inputs["big_inlet_speed"] = round(randfp(0.2, 0.3, precision) / (2 * math.pi * inputs["r4"]), precision)
  inputs["small_inlet_speed"] = round(randfp(0.065, 0.07, precision) * inputs["big_inlet_speed"], precision) * -1.0 # Negative z direction
  inputs["kinematic_viscosity"] = round(0.00001789 * 287.057 * inputs["temperature"] / p_atm, 14)
monte_carlo() # Refresh to initialy add physical parameters to inputs dictionary

def save_inputs(filename):
  """
  Save all input parameters to a JSON file.
  """
  with open(filename, 'w') as f:
    json.dump(inputs, f, indent=2)

# MESH GENERATION

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
  #gmsh.model.getEntities(dim=-1) # Print all entities
  # Merge all volumes into a single one
  volumes = [volume1, volume2, volume3, volume4]
  volume_tuples = [(3, v) for v in volumes]
  fused_volume, _ = gmsh.model.occ.fuse(volume_tuples, volume_tuples)
  gmsh.model.occ.synchronize() # Sync changes
  # Physical boundaries
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
  # Mesh size fields for small inlet
  gmsh.model.mesh.field.add("Distance", 1)
  gmsh.model.mesh.field.setNumbers(1, "FacesList", [4, 6])
  #gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", inputs["points_per_curve"])
  gmsh.model.mesh.field.add("Threshold", 2)
  gmsh.model.mesh.field.setNumber(2, "InField", 1)
  gmsh.model.mesh.field.setNumber(2, "SizeMin", inputs["small_inlet_minimum_mesh_size"]) # Minimum element size in the small inlet
  gmsh.model.mesh.field.setNumber(2, "SizeMax", inputs["small_inlet_maximum_mesh_size"]) # Maximum element size elsewhere
  gmsh.model.mesh.field.setNumber(2, "DistMin", inputs["small_inlet_minimum_distance"])
  gmsh.model.mesh.field.setNumber(2, "DistMax", inputs["small_inlet_maximum_distance"])
  gmsh.model.mesh.field.setAsBackgroundMesh(2)
  # Set options
  gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
  gmsh.option.setNumber("General.Terminal", 1) # Enable terminal output
  gmsh.option.setNumber("Mesh.Algorithm", 6) # Front-Delaunay default algorithm
  gmsh.option.setNumber("General.NumThreads", os.cpu_count()) # Parallel processing with all threads
  if inputs["mesh_type"]:  # quad mesh specific options
    gmsh.option.setNumber("Mesh.RecombineAll", 1) # Triangular recombination into quadrilaterals
    surfaces = gmsh.model.getEntities(dim=2)
    for surface in surfaces:
      gmsh.model.mesh.setRecombine(2, surface[1])
      gmsh.model.mesh.setSmoothing(2, surface[1], 10) # Smoothing to improve element quality
  gmsh.model.mesh.generate(3) # Generate 3D mesh
  #gmsh.write("fluid.msh") # Write mesh to file
  mesh_file = f"fluid_smax{inputs['small_inlet_maximum_mesh_size']}_smin{inputs['small_inlet_minimum_mesh_size']}_pc{inputs['points_per_curve']}_dmax{inputs['small_inlet_maximum_distance']}_dmin{inputs['small_inlet_minimum_distance']}.msh"
  gmsh.write(mesh_file) # Write mesh to file
  # gmsh.fltk.run() # Visualize
  gmsh.finalize()
  return mesh_file

# CREATE SIMULATION CASE FROM SCRATCH

# Configuration arrays for case set up
header = ["standalone", """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  11
     \\\\/     M anipulation  |
\*---------------------------------------------------------------------------*/\n"""]
line_mid = ["standalone", "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"]
line_end = ["standalone", "// ************************************************************************* //"]
new_line = ["standalone", ""]

p = [header,
    ["block", "FoamFile", [
      ["parameter", "format", "ascii"],
      ["parameter", "class", "volScalarField"],
      ["parameter", "object", "p"]
    ]],
    line_mid,
    ["parameter", "dimensions", "[0 2 -2 0 0 0 0]"],
    ["parameter", "internalField", "uniform 0"],
    ["block", "boundaryField", [
      ["block", "big-inlet", [
        ["parameter", "type", "zeroGradient"]
      ]],
      ["block", "big-outlet", [
        ["parameter", "type", "fixedValue"],
        ["parameter", "value", "uniform 0"]
      ]],
      ["block", "small-inlet", [
        ["parameter", "type", "zeroGradient"]
      ]],
      ["block", "walls", [
        ["parameter", "type", "zeroGradient"]
      ]]
    ]],
    line_end]

U = [header,
    ["block", "FoamFile", [
      ["parameter", "format", "ascii"],
      ["parameter", "class", "volVectorField"],
      ["parameter", "object", "U"]
    ]],
    line_mid,
    ["parameter", "dimensions", "[0 1 -1 0 0 0 0]"],
    ["parameter", "internalField", "uniform (0 0 0)"],
    ["block", "boundaryField", [
      ["block", "big-inlet", [
        ["parameter", "type", "fixedValue"],
        ["parameter", "value", f"uniform ({inputs['big_inlet_speed']} 0 0)"]
      ]],
      ["block", "big-outlet", [
        ["parameter", "type", "zeroGradient"]
      ]],
      ["block", "small-inlet", [
        ["parameter", "type", "fixedValue"],
        ["parameter", "value", f"uniform (0 0 {inputs['small_inlet_speed']})"]
      ]],
      ["block", "walls", [
        ["parameter", "type", "noSlip"]
      ]]
    ]],
    line_end]

physicalProperties = [header,
    ["block", "FoamFile", [
      ["parameter", "format", "ascii"],
      ["parameter", "class", "dictionary"],
      ["parameter", "location", "\"constant\""],
      ["parameter", "object", "physicalProperties"]
    ]],
    line_mid,
    ["parameter", "nu", f"[0 2 -1 0 0 0 0] {inputs['kinematic_viscosity']}"],
    new_line,
    line_end]

controlDict = [header,
    ["block", "FoamFile", [
      ["parameter", "format", "ascii"],
      ["parameter", "class", "dictionary"],
      ["parameter", "location", "\"system\""],
      ["parameter", "object", "controlDict"]
    ]],
    line_mid,
    ["parameter", "application", "icoFoam"],
    ["parameter", "startFrom", "latestTime"],
    ["parameter", "startTime", f"{inputs['starttime']}"],
    ["parameter", "stopAt", "endTime"],
    ["parameter", "endTime", f"{inputs['endtime']}"],
    ["parameter", "deltaT", f"{inputs['deltat']}"],
    ["parameter", "writeControl", f"{inputs['writeControl']}"],
    ["parameter", "writeInterval", f"{inputs['writeinterval']}"],
    ["parameter", "purgeWrite", "0"],
    ["parameter", "writeFormat", "ascii"],
    ["parameter", "writePrecision", "6"],
    ["parameter", "writeCompression", "off"],
    ["parameter", "timeFormat", "general"],
    ["parameter", "timePrecision", "6"],
    ["parameter", "runTimeModifiable", "true"],
    line_end]

foamDataToFluentDict = [header,
    ["block", "FoamFile", [
      ["parameter", "format", "ascii"],
      ["parameter", "class", "dictionary"],
      ["parameter", "location", "\"system\""],
      ["parameter", "object", "foamDataToFluentDict"]
    ]],
    line_mid,
    ["parameter", "p", "1"],
    ["parameter", "U", "2"],
    ["parameter", "T", "3"],
    ["parameter", "h", "4"],
    ["parameter", "k", "5"],
    ["parameter", "epsilon", "6"],
    ["parameter", "alpha1", "150"],
    new_line,
    line_end]

fvSchemes = [header,
    ["block", "FoamFile", [
      ["parameter", "format", "ascii"],
      ["parameter", "class", "dictionary"],
      ["parameter", "location", "\"system\""],
      ["parameter", "object", "fvSchemes"]
    ]],
    line_mid,
    ["block", "ddtSchemes",[
      ["parameter", "default", "Euler"]
    ]],
    ["block", "gradSchemes",[
      ["parameter", "default", "Gauss linear"]
    ]],
    ["block", "divSchemes", [
      ["parameter", "default", "none"],
      ["parameter", "div(phi,U)", "Gauss limitedLinearV 1"]
    ]],
    ["block", "laplacianSchemes",[
      ["parameter", "default", "Gauss linear corrected"]
    ]],
    ["block", "interpolationSchemes",[
      ["parameter", "default", "linear"]
    ]],
    ["block", "snGradSchemes",[
      ["parameter", "default", "corrected"]
    ]],
    new_line,
    line_end]

fvSolution = [header,
    ["block", "FoamFile", [
      ["parameter", "format", "ascii"],
      ["parameter", "class", "dictionary"],
      ["parameter", "location", "\"system\""],
      ["parameter", "object", "fvSolution"]
    ]],
    line_mid,
    ["block", "solvers", [
      ["block", "p", [
        ["parameter", "solver", "PCG"],
        ["parameter", "preconditioner", "DIC"],
        ["parameter", "tolerance", "1e-06"],
        ["parameter", "relTol", "0.05"]
      ]],
      ["block", "pFinal", [
        ["parameter", "$p", ""],
        ["parameter", "relTol", "0"]
      ]],
      ["block", "U", [
        ["parameter", "solver", "smoothSolver"],
        ["parameter", "smoother", "symGaussSeidel"],
        ["parameter", "tolerance", "1e-05"],
        ["parameter", "relTol", "0"]
      ]],
    ]],
    ["block", "PISO", [
      ["parameter", "nCorrectors", "2"],
      ["parameter", "nNonOrthogonalCorrectors", "2"]
    ]],
    new_line,
    line_end]

def dir(base_dir, prefix='s'):
    """
    Create a new directory inside "base_dir" named "prefix + number", if "prefix + number - 1" already exists.
    """
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

def create_file_structure(file_paths):
  """
  Creates files and directories based on the input array.
  args:
    file_paths (list of lists): [[(str) type, (str) path],... ]
      type can be: "f" for file and "d" for directory
  """
  for item in file_paths:
    if len(item) != 2:
      raise ValueError(f"Invalid input '{item}'. Each inner list should have exactly 2 elements.")
      continue
    file_type, path = item
    if file_type == "f": # Create a file
      try:
        with open(path, "w"):
          pass  # Create an empty file
        print(f"Created file:  {path}")
      except Exception as e:
        print(f"Error creating file {path}: {e}")
    elif file_type == "d": # Create a directory
      try:
        os.makedirs(path, exist_ok=True)
        print(f"Created directory:  {path}")
      except Exception as e:
        print(f"Error creating directory {path}: {e}")
    else:
      print(f"Invalid filetype: {file_type}. Use 'f' for file or 'd' for directory.")

def param2text(param, indent):
  """
  Parser function to generate OpenFOAM's dictionaries and configuration files from configuration arrays.
  args:
    param (list): [ (str) type, (str) name, (str, num or list) value]
      type can be: "standalone", "parameter" and "block"
      if type is standalone: ["standalone", value]
      if type is block, value is a list of types (elements)
    indent (integer >= 0)
  return:
    text: generated string
  """
  text = ""
  tabs = indent * "\t"
  match param[0]: # Check parameter type
    case "standalone":
      if len(param) != 2:
        raise ValueError("Parameter type \"standalone\" expects 2 elements. Received: " + str(param))
      text += param[1]
    case "parameter":
      if len(param) != 3:
        raise ValueError("Parameter type \"parameter\" expects 3 elements. Received: " + str(param))
      if param[2] == "": # if paramater has no value (only name)
        text += tabs + param[1] + ";\n"
      else: # if parameter is a name-value data pair
        text += tabs + param[1] + "\t" + param[2] + ";\n"
    case "block":
      if len(param) != 3:
        raise ValueError("Parameter type \"block\" expects 3 elements. Received: " + str(param))
      text += tabs + param[1] + "\n" + tabs + "{\n" # Block name {
      for prmt in param[2]:
        text += param2text(prmt, indent+1) # Block contents
      text += tabs + "}\n" # } End block
    case _: # Unkown type
      raise ValueError("Unknown parameter type: " + str(param[0]))
  return text

def populate_file(file_name, file_path):
  """
  Given a file name and its path, generate content (a configuration) for every simulation file.
  """
  match file_name:
    case "p":
      content = p
    case "U":
      content = U
    case "physicalProperties":
      content = physicalProperties
    case "controlDict":
      content = controlDict
    case "foamDataToFluentDict":
      content = foamDataToFluentDict
    case "fvSchemes":
      content = fvSchemes
    case "fvSolution":
      content = fvSolution
    case _: # Unkown file name
      raise ValueError(f"The given file name '{file_name}' is unknown so no configuration can be generated for it.")
  try:
    with open(file_path, 'a') as config_file:
      for param in content:
        config_file.write(param2text(param,0))
    print(f"Configuration written to {file_path}")
  except FileNotFoundError:
    print(f"Unreachable path '{file_path}'.")
    print(param2text(param, 0))
    return 1
  return 0

def build_case(base_path = "."):
  """
  Creates the directories and files necessary for a simulation and populates them.
  """
  file_structure = [
    ["d", base_path + "/0"],
    ["f", base_path + "/0/p"],
    ["f", base_path + "/0/U"],
    ["d", base_path + "/constant"],
    ["f", base_path + "/constant/physicalProperties"],
    ["d", base_path + "/system"],
    ["f", base_path + "/system/controlDict"],
    ["f", base_path + "/system/foamDataToFluentDict"],
    ["f", base_path + "/system/fvSchemes"],
    ["f", base_path + "/system/fvSolution"]
  ]
  create_file_structure(file_structure)
  for file in file_structure:
    if file[0] == "f": # skip directories
      file_name = file[1].split("/")[-1] # extract file name
      populate_file(file_name, file[1])
  print()

# MODIFY EXISTING CASE

def cpdir(base_dir, source_dir, prefix='s'):
  """
  Function to recursive copy directories and their contents.
  Used to copy template cases.
  """
  counter = 0
  while True:
    dir_name = f"{prefix}{counter}"
    dir_path = os.path.join(base_dir, dir_name)
    if not os.path.exists(dir_path):
      print(f"Copied directory:  '{source_dir}'->'{dir_name}'")
      shutil.copytree(source_dir, dir_path)
      return dir_name
    else:
        counter += 1

def mod_case(path):
  """
  Modify existing case files using PyFoam.RunDictionary.ParsedParameterFile function.
  """
  # p
  p_file = ParsedParameterFile(f'{path}/0/p')
  #p_file["boundaryField"][boundary_conditions[0]]["type"] = "zeroGradient" 
  #p_file["boundaryField"][boundary_conditions[1]]["type"] = "fixedValue"
  p_file["boundaryField"][inputs["boundary_conditions"][1]]["value"] = "uniform 0"
  #p_file["boundaryField"][boundary_conditions[2]]["type"] = "zeroGradient"
  #p_file["boundaryField"][boundary_conditions[3]]["type"] = "zeroGradient"
  p_file.writeFile()
  print(f"Modified pressure parameters of file '{path}/0/p'.")
  # U
  U_file = ParsedParameterFile(f'{path}/0/U')
  #U_file["boundaryField"][boundary_conditions[0]]["type"] = "fixedValue" 
  U_file["boundaryField"][inputs["boundary_conditions"][0]]["value"] = f"uniform ({inputs['big_inlet_speed']} 0 0)"
  #U_file["boundaryField"][boundary_conditions[1]]["type"] = "zeroGradient"
  #U_file["boundaryField"][boundary_conditions[2]]["type"] = "fixedValue"
  U_file["boundaryField"][inputs["boundary_conditions"][2]]["value"] = f"uniform (0 0 {inputs['small_inlet_speed']})"
  #U_file["boundaryField"][boundary_conditions[3]]["type"] = "noSlip"
  U_file.writeFile()
  print(f"Modified velocity parameters of file '{path}/0/U'.")
  # physicalProperties
  physical_properties = ParsedParameterFile(f'{path}/constant/physicalProperties')
  physical_properties["nu"] = f"[0 2 -1 0 0 0 0] {inputs['kinematic_viscosity']}"
  # controlDict
  #control_dict = ParsedParameterFile(f'{path}/system/controlDict')
  #control_dict["application"] = "icoFoam"
  #control_dict["startTime"] = starttime
  #control_dict["endTime"] = endtime
  #control_dict["deltaT"] = deltat
  #control_dict["writeInterval"] = writeinterval
  #control_dict.writeFile()
  #print(f"Modified controlDict parameters of file '{path}/system/controlDict'.")
  print()

# RESULT EXTRACTION AND PROCESSING

def csv_template(csv_path, header):
  """
  Given a path and a header, creates a CSV file in path with the corresponding header
  """
  with open(csv_path, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header)

def vtk_results(sn, vtk_dir, csv_path, x_range, y_range, z_range, area):
  """
  Receives:
    * a simulation number (sn),
    * directory containing VTK files (vtk_dir),
    * a target CSV file to store the resulst (csv_path),
    * and optionally tupples contianing the x, y, and z ranges of the geometrical space to limit the data to.
  Returns nothing.
  Appends the data extracted from a single simulation's VTK files to the global simulation result CSV, containing all simulation results.
  """
  vtk_files = [f for f in os.listdir(vtk_dir) if f.endswith('.vtk')]
  if not vtk_files:
    print(f"No VTK files found in the directory: {vtk_dir}")
    return
  vtk_files.sort(key=lambda f: int(f.split('_')[-1].replace('.vtk', '')))
  for vtk_file in vtk_files:
    vtk_path = os.path.join(vtk_dir, vtk_file)
    print(f"Reading VTK file '{vtk_path}'")
    vtk_lector = vtk.vtkUnstructuredGridReader()
    vtk_lector.SetFileName(vtk_path)
    vtk_lector.Update()
    vtk_data = vtk_lector.GetOutput()
    if not vtk_data:
      print("Failed to read VTK file.")
      return
    if not vtk_data.IsA('vtkUnstructuredGrid'):
      print(f"Unsupported VTK data type in file: {vtk_file}")
    else:
      points = vtk_data.GetPoints()
      if not points:
        print("No points found in the VTK data.")
        return
      num_points = points.GetNumberOfPoints()
      points_array = numpy.zeros((num_points, 3))
      for i in range(num_points):
        points_array[i] = points.GetPoint(i)
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
      time_step = vtk_file.split('_')[-1].replace('.vtk', '')
      with open(csv_path, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for i in range(num_points):
          x, y, z = points_array[i]
          # Filter points geometrically wise
          if x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1] and z_range[0] <= z <= z_range[1]:
            row = [
              sn,
              time_step,
              x, y, z,
              pressure_values[i],
              velocity_values[i][0], velocity_values[i][1], velocity_values[i][2],
              area
            ]
            csvwriter.writerow(row)
  print(f"\nResults for simulation {sn} written to '{csv_path}'")

# Material parameters
# Prices of Bright Drawn Bar 304 Austenitic Stainless Steel in Europe (frebruary 2024)
material_cost = 2.999 # EUR (€) / kg
material_density = 7930.0 # kg / m³
# cost (€/mm) = area (m²) * material_cost (€/kg) * material_density (kg/m³) * unit_conversion_factor (10⁶mm²/1m² * 1m³/10⁹mm³ = 1m/10³mm = 0.001 m/mm)
  
def collect_results(input_csv, output_csv):
  """
  Collect results to a CSV file containing each simulation point (average pressure, cost per thickness).
  Calulate average pressures (Pa) and cost per component thickness (€/mm)
  """
  # Load the data from the CSV file
  data = pd.read_csv(input_csv)
  # Sort the dataframe to ensure the latest time comes last in each group
  data.sort_values(by=['sn', 'time'], inplace=True)
  # Group by 'sn' to process each simulation
  grouped = data.groupby('sn')
  # Initialize an empty DataFrame for results with the appropriate columns
  result = pd.DataFrame(columns=['sn', 'cost', 'avg_p'])
  result = result.astype({'sn': 'int', 'cost': 'float', 'avg_p': 'float'})
  for name, group in grouped:
    # Get the last time step's data
    last_time_step = group[group['time'] == group['time'].max()]
    # Calculate the average pressure of the last time step
    avg_p = last_time_step['p'].mean()
    # As area is constant for each simulation take the first one
    area = group['area'].iloc[0]
    # Compute the value using material cost and density
    cost = material_cost * material_density * 0.001 * area
    # Create a new DataFrame for this iteration's result
    new_result = pd.DataFrame({'sn': [name], 'cost': [cost], 'avg_p': [avg_p]})
    new_result = new_result.astype({'sn': 'int', 'cost': 'float', 'avg_p': 'float'})
    # Ensure no all-NA or empty columns are present before concatenation
    if not new_result.empty:
      result = pd.concat([result, new_result], ignore_index=True)
  # Write the result to a new CSV file
  result.to_csv(output_csv, index=False)

# PARETO

def not_dominated(point, others):
  """
  Check if a data point is dominated by others or a Pareto optimal point.
  """
  # A point is dominated if there is any other point with lower 'avg_p' and 'cost'
  for other in others:
    if other['avg_p'] <= point['avg_p'] and other['cost'] <= point['cost'] and (other['avg_p'] < point['avg_p'] or other['cost'] < point['cost']):
      return False
  return True

def find_pareto_points(input_csv, output_csv):
  """
  Given an input CSV, write Pareto optimal points to output CSV.
  """
  # Load CSV data
  data = pd.read_csv(input_csv)
  pareto_points = []
  # Convert DataFrame to list of dicts for easier manipulation
  points = data.to_dict('records')
  # Check for Pareto optimality
  for point in points:
    if not_dominated(point, points):
      pareto_points.append(point)
  # Convert Pareto points back to DataFrame for easier viewing/manipulation
  pareto_df = pd.DataFrame(pareto_points)
  # Save the DataFrame of Pareto points to a new CSV file
  pareto_df.to_csv(output_csv, index=False)
  return pareto_df, data

# PLOT

def pareto(input_csv, output_csv):
  """
  Compute Pareto optimal points from input CSV,
  write points to output CSV,
  and generate the three following plots:
    * all points, Pareto points and frontier
    * all points and Pareto points highlighted
    * only Pareto points
    * zoomed plot of Pareto optimal points with highlighted frontier
  """
  pareto_df, all_data = find_pareto_points(input_csv, output_csv)
  # Plot 1: All points, Pareto points and frontier
  plt.figure(figsize=(10, 6))
  plt.scatter(all_data['cost'], all_data['avg_p'], color='blue', marker='o', label='All Points')
  plt.scatter(pareto_df['cost'], pareto_df['avg_p'], color='red', marker='o', label='Pareto Points')
  pareto_sorted = pareto_df.sort_values(by='cost') # Sort points to draw Pareto frontier
  plt.plot(pareto_sorted['cost'], pareto_sorted['avg_p'], color='green', label='Pareto Frontier')
  plt.title('All Points vs. Pareto Optimal Points with Frontier')
  plt.xlabel('Cost [€/mm]')
  plt.ylabel('Average Pressure [Pa]')
  plt.legend()
  plt.savefig('plot_all_pareto_frontier.pdf')
  plt.show()
  # Plot 2: All points and Pareto points highlighted
  plt.figure(figsize=(10, 6))
  plt.scatter(all_data['cost'], all_data['avg_p'], color='blue', marker='o', label='All Points')
  plt.scatter(pareto_df['cost'], pareto_df['avg_p'], color='red', marker='o', label='Pareto Points')
  plt.title('All Points vs. Pareto Optimal Points')
  plt.xlabel('Cost [€/mm]')
  plt.ylabel('Average Pressure [Pa]')
  plt.legend()
  plt.savefig('plot_pareto_all.pdf')  # Save the plot to a file
  plt.close()  # Close the plot to free up memory
  # Plot 3: Only Pareto points
  plt.figure(figsize=(10, 6))
  plt.scatter(pareto_df['cost'], pareto_df['avg_p'], color='red', marker='o', label='Pareto Points')
  plt.title('Pareto Optimal Points')
  plt.xlabel('Cost [€/mm]')
  plt.ylabel('Average Pressure [Pa]')
  plt.legend()
  plt.savefig('plot_pareto.pdf')  # Save the plot to a file
  plt.close()  # Close the plot to free up memory
  # Plot 4: Zoomed-in plot of Pareto optimal points with highlighted frontier
  plt.figure(figsize=(10, 6))
  plt.scatter(all_data['cost'], all_data['avg_p'], color='blue', marker='o', alpha=0.3, label='All Points')
  plt.scatter(pareto_df['cost'], pareto_df['avg_p'], color='red', marker='o', label='Pareto Points')
  plt.plot(pareto_sorted['cost'], pareto_sorted['avg_p'], color='green', label='Pareto Frontier')
  # Set axis limits to zoom in on Pareto points
  padding = 0.1  # Add a 10% margin around Pareto points
  x_min, x_max = pareto_df['cost'].min(), pareto_df['cost'].max()
  y_min, y_max = pareto_df['avg_p'].min(), pareto_df['avg_p'].max()
  plt.xlim(x_min - padding * (x_max - x_min), x_max + padding * (x_max - x_min))
  plt.ylim(y_min - padding * (y_max - y_min), y_max + padding * (y_max - y_min))
  plt.title('Pareto Optimal Points with Frontier')
  plt.xlabel('Cost [€/mm]')
  plt.ylabel('Average Pressure [Pa]')
  plt.legend()
  plt.savefig('plot_all_zoomed_pareto_frontier.pdf')
  plt.close()

# STATISTICS

def plots_cost_pressure(input_csv):
  # Initialize lists to hold data
  sn_list = []
  costs = []
  avg_pressures = []
  # Read data from the CSV file
  with open(input_csv, mode='r', newline='') as file:
    reader = csv.DictReader(file)
    for row in reader:
      sn = int(row['sn'])
      cost = float(row['cost'])
      avg_p = float(row['avg_p'])
      # Append data to lists
      sn_list.append(sn)
      costs.append(cost)
      avg_pressures.append(avg_p)
  # Plot 1: Cost vs Simulation Number
  plt.figure(figsize=(10, 6))
  plt.scatter(sn_list, costs, color='blue', marker='o', label='Cost vs Simulation Number')
  plt.title('Cost vs Simulation Number')
  plt.xlabel('Simulation Number')
  plt.ylabel('Cost [€/mm]')
  plt.grid(True)
  plt.legend()
  plt.savefig("plot_sn_cost.pdf")
  plt.close()
  # Plot 2: Average Pressure vs Simulation Number
  plt.figure(figsize=(10, 6))
  plt.scatter(sn_list, avg_pressures, color='blue', marker='o', label='Avg Pressure vs Simulation Number')
  plt.title('Average Pressure vs Simulation Number')
  plt.xlabel('Simulation Number')
  plt.ylabel('Average Pressure [Pa]')
  plt.grid(True)
  plt.legend()
  plt.savefig("plot_sn_pressure.pdf")
  plt.close()
  # Plot 3: Average Pressure vs Cost
  plt.figure(figsize=(10, 6))
  plt.scatter(costs, avg_pressures, color='blue', marker='o', label='Pressure vs Cost')
  plt.title('Average Pressure vs Cost')
  plt.xlabel('Cost [€/mm]')
  plt.ylabel('Average Pressure [Pa]')
  plt.grid(True)
  plt.legend()
  plt.savefig("plot_cost_pressure.pdf")
  plt.close()

def compute_statistics(csv_file):
  """
  Given a CSV file containing the simulation number, cost and average inlet pressure,
  compute statistical metrics and store them on a JSON file.
  """
  # Initialize lists to store data
  sn_list = []
  cost_list = []
  avg_p_list = []
  # Load CSV data
  with open(csv_file, 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
      sn_list.append(int(row['sn']))
      cost_list.append(float(row['cost']))
      avg_p_list.append(float(row['avg_p']))
  # Compute statistical metrics
  statistics_data = {
    "number_of_simulations": len(sn_list),
    "cost": {
      "max": max(cost_list),
      "min": min(cost_list),
      "average": sum(cost_list) / len(cost_list),
      "median": statistics.median(cost_list),
      #"mode": statistics.mode(area_list) if len(set(area_list)) != len(area_list) else "No mode"
      "std_dev": statistics.stdev(cost_list),
      "variance": statistics.variance(cost_list),
    },
    "avg_p": {
      "max": max(avg_p_list),
      "min": min(avg_p_list),
      "average": sum(avg_p_list) / len(avg_p_list),
      "median": statistics.median(avg_p_list),
      #"mode": statistics.mode(avg_p_list) if len(set(avg_p_list)) != len(avg_p_list) else "No mode"
      "std_dev": statistics.stdev(avg_p_list),
      "variance": statistics.variance(avg_p_list),
    }
  }
  # Place timestamp on filename to prevent overwriting
  json_filename = f"statistics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
  # Write the data to a JSON file
  with open(json_filename, 'w') as json_file:
    json.dump(statistics_data, json_file, indent=2)
  return json_filename

