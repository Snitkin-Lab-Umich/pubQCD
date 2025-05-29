__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
import json
from snakemake.shell import shell

# Retrieve params and paths
size = snakemake.params.get("size", "")
coverage_output = snakemake.output.coverage

def run_fastq_scan():
    shell("zcat {snakemake.input.r1} {snakemake.input.r2} | fastq-scan -g {size} > {coverage_output}")

try:
    run_fastq_scan()
except Exception as e1:
    print(f"First attempt to run fastq-scan failed: {e1}")
    try:
        print("Retrying fastq-scan...")
        run_fastq_scan()
    except Exception as e2:
        print(f"Second attempt to run fastq-scan failed: {e2}")
        print("Writing dummy JSON output.")
        dummy_data = {
            "qc_stats": {
                "total_bp": 0,
                "coverage": 0,
                "read_total": 0,
                "read_mean": 0,
            }
        }
        with open(coverage_output, "w") as f:
            json.dump(dummy_data, f, indent=4)