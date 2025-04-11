__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from os import path
import os
import ruamel.yaml

# Use snakemake params to get base_dir and prefix
base_dir = snakemake.params.get("base_dir", "")
prefix = snakemake.params.get("prefix", "")

# Define the directories based on provided paths
pipeline_dir = base_dir if base_dir else os.path.dirname(os.path.abspath(__file__))
QCD_dir = os.path.dirname(pipeline_dir)
# Paths to configuration files
config_file = os.path.join(QCD_dir, "config", "config_pass_fail.yaml")
samples_passed_coverage_file = os.path.join(QCD_dir, "results", prefix, "sample_files", "samples_passed_coverage.csv")
update_path_success_file = os.path.join(QCD_dir, "results", prefix, "sample_files", "updated_path.done")

# Ensure the file path is normalized
config_file = path.normpath(config_file)
samples_passed_coverage_file = path.normpath(samples_passed_coverage_file)

# Load the existing YAML content
yaml = ruamel.yaml.YAML(typ="rt")
yaml.width = 100

# Debug logging
print(f"Loading config from: {config_file}")
print(f"Updating 'samples_passed_coverage' path to: {samples_passed_coverage_file}")

with open(config_file, "r") as file:
    config = yaml.load(file)

# Update if the key is present
if "samples_passed_coverage" in config:
    config["samples_passed_coverage"] = samples_passed_coverage_file

    # Write the updated config back to the file
    with open(config_file, "w") as file:
        yaml.dump(config, file)

    # Create the success marker file
    with open(update_path_success_file, "w") as done_file:
        done_file.write(f"Path:{samples_passed_coverage_file} updated successfully.\n")
    
    print(f"Path updated successfully here {update_path_success_file}")

else:
    raise KeyError("'samples_passed_coverage' key not found in the config file.")