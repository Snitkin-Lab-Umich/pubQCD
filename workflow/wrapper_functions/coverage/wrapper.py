__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from os import path
import os
from snakemake.shell import shell
import pandas as pd
import json

# import ruamel.yaml

size = snakemake.params.get("size", "")
outdir = snakemake.params.get("outdir", "")
# r1_files = snakemake.input.r1  
# r2_files = snakemake.input.r2 

passed_csv=os.path.join(outdir,"sample_files", "samples_passed_coverage.csv")

shell("zcat {snakemake.input.r1} {snakemake.input.r2} | fastq-scan -g {size} > {snakemake.output.coverage}")

with open({snakemake.output.coverage}) as f:
        data = json.load(f)

sample_name = os.path.basename(coverage_json).replace("_coverage.json", "")
failed_csv = os.path.join(outdir, "sample_files/samples_failed_coverage_summary.csv")

coverage = data.get('qc_stats', {}).get('coverage', 0)
total_reads = data.get('qc_stats', {}).get('read_total', 0)
total_bp = data.get('qc_stats', {}).get('total_bp', 0)
mean_read_length = data.get('qc_stats', {}).get('read_mean', 0)

# Ensure passed_csv exists before appending
if not os.path.exists((passed_csv)):
    pd.DataFrame(columns=["sample_id", "illumina_r1"]).to_csv(passed_csv, index=False)

# Read existing passed samples to avoid duplicates
existing_passed_samples = set(pd.read_csv(passed_csv)["sample_id"])

if coverage > 20 and sample_name not in existing_passed_samples:
    passed_df = pd.DataFrame([[sample_name, f"{sample_name}_R1.fastq.gz"]],
                                columns=["sample_id", "illumina_r1"])
    passed_df.to_csv((passed_csv), mode='a', index=False, header=False)

    # with open(updated_samples_file, "a") as done_file:
    #     done_file.write(f"Sample: {sample_name} updated successfully.\n")

elif coverage <= 20:
    failed_df = pd.DataFrame([[sample_name, total_reads, total_bp, mean_read_length, coverage] + ["NA"] * 15 + ["FAIL"] + ["NA"] * 5],
                                columns=[
                                    "Sample", "Total_reads", "Total_bp", "MeanReadLength", "Coverage", "Scheme",
                                    "ST", "After_trim_per_base_sequence_content", "After_trim_overrepresented_sequences",
                                    "After_trim_%GC", "After_trim_Total Bases", "After_trim_Total Sequences",
                                    "After_trim_median_sequence_length", "After_trim_avg_sequence_length",
                                    "After_trim_total_deduplicated_percentage", "After_trim_Sequence length",
                                    "After_trim_adapter_content", "N50", "Total length", "Total # of contigs",
                                    "QC Check", "ANI", "Align_fraction_ref", "Align_fraction_query",
                                    "Ref_name", "Species"
                                ])
    failed_df.to_csv(failed_csv, mode='a', index=False, header=not os.path.exists(failed_csv))


#shell("samtools view -Sb {snakemake.output.clipped_sam_out} > {snakemake.output.bam_out} && samtools sort {snakemake.output.bam_out} -m 500M -@ 0 -o {snakemake.output.sorted_bam_out} -T {outdir_temp} &> {snakemake.log.post_align_sam_to_bam_log}")

# Use snakemake params to get base_dir and prefix
# prefix = snakemake.params.get("prefix", "")

# # Define the directories based on provided paths
# pipeline_dir = base_dir if base_dir else os.path.dirname(os.path.abspath(__file__))
# QCD_dir = os.path.dirname(pipeline_dir)
# # Paths to configuration files
# config_file = os.path.join(QCD_dir, "config", "config_pass_fail.yaml")
# samples_passed_coverage_file = os.path.join(QCD_dir, "results", prefix, "sample_files", "samples_passed_coverage.csv")
# update_path_success_file = os.path.join(QCD_dir, "results", prefix, "sample_files", "updated_path.done")

# # Ensure the file path is normalized
# config_file = path.normpath(config_file)
# samples_passed_coverage_file = path.normpath(samples_passed_coverage_file)

# # Load the existing YAML content
# yaml = ruamel.yaml.YAML(typ="rt")
# yaml.width = 100

# # Debug logging
# print(f"Loading config from: {config_file}")
# print(f"Updating 'samples_passed_coverage' path to: {samples_passed_coverage_file}")

# with open(config_file, "r") as file:
#     config = yaml.load(file)

# # Update if the key is present
# if "samples_passed_coverage" in config:
#     config["samples_passed_coverage"] = samples_passed_coverage_file

#     # Write the updated config back to the file
#     with open(config_file, "w") as file:
#         yaml.dump(config, file)

#     # Create the success marker file
#     with open(update_path_success_file, "w") as done_file:
#         done_file.write(f"Path:{samples_passed_coverage_file} updated successfully.\n")
    
#     print(f"Path updated successfully here {update_path_success_file}")

# else:
#     raise KeyError("'samples_passed_coverage' key not found in the config file.")