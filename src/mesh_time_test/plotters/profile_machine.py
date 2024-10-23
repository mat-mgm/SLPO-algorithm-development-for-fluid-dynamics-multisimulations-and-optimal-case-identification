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

import sys
if len(sys.argv) > 1:
  print("The first argument must be a path to a JSON file with mesh size and timming data. All other arguments provided will be ignored.")
  profile_machine(sys.argv[1])
else:
  print("No argument was provided.")

