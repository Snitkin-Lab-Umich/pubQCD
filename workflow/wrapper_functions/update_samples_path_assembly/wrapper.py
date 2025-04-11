__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2025, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from os import path
import os
import sys
import shutil
import ruamel.yaml

# Use snakemake params to get base_dir and prefix
base_dir = snakemake.params.get("base_dir", "")
prefix = snakemake.params.get("prefix", "")

# Define the directories based on provided paths
pipeline_dir = base_dir if base_dir else os.path.dirname(os.path.abspath(__file__))
QCD_dir = os.path.dirname(pipeline_dir)
# Paths to configuration files
config_file = os.path.join(QCD_dir, "config", "config_pass_fail.yaml")

samples_passed_assembly_step = os.path.join(QCD_dir, "results", prefix, "sample_files", "samples_passed_assembly_step.csv") # Samples passed coverage file from assembly step
samples_passed_assembly_file = os.path.join(QCD_dir, "results", prefix, "sample_files", "samples_passed_assembly.csv")

# Ensure source file exists before copying
if os.path.exists(samples_passed_assembly_step):
    shutil.copyfile(samples_passed_assembly_step, samples_passed_assembly_file)
    print(f"Copied: {samples_passed_assembly_step} to {samples_passed_assembly_file}")
else:
    print(f"Error: Source file {samples_passed_assembly_step} does not exist!")
    sys.exit(1)

update_path_success_file = os.path.join(QCD_dir, "results", prefix, "sample_files", "updated_path_assembly.done")

# Ensure the file path is normalized
config_file = path.normpath(config_file)
samples_passed_assembly_file = path.normpath(samples_passed_assembly_file)

# Load the existing YAML content
yaml = ruamel.yaml.YAML(typ="rt")
yaml.width = 100

# Debug logging
print(f"Loading config from: {config_file}")
print(f"Updating 'samples_passed_assembly' path to: {samples_passed_assembly_file}")

with open(config_file, "r") as file:
    config = yaml.load(file)

# Update if the key is present
if "samples_passed_assembly" in config:
    config["samples_passed_assembly"] = samples_passed_assembly_file

    # Write the updated config back to the file
    with open(config_file, "w") as file:
        yaml.dump(config, file)

    # Create the success marker file
    with open(update_path_success_file, "w") as done_file:
        done_file.write(f"Path:{samples_passed_assembly_file} updated successfully.\n")
    
    print(f"Path updated successfully here {update_path_success_file}")

else:
    raise KeyError("'samples_passed_assembly' key not found in the config file.")