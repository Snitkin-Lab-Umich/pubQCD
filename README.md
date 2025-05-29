# pubQCD - Public datasets Quality Control and Contamination Detection workflow

pubQCD is a quality control pipeline for datasets downloaded from public database implemented in Snakemake that takes either raw fastqs or assemblies as input.

If you downloaded assemblies, click [here](README_assemblies.md), otherwise, continue reading the instructions on how to process raw reads. **You do not need to download the fastqs as the pipeline does it for you.**

### Summary

As part of the SOP in the [Snitkin lab](https://thesnitkinlab.com/index.php), this pipeline should be run on raw sequencing data (separate pipeline for assemblies coming soon!) as soon as the data is downloaded from a public database. 

In short, it performs the following steps:

* [Fastqc](https://github.com/s-andrews/FastQC) is used to generate HTML reports to asses quality of sequencing reads before and after trimming reads. 
* Trims and filters low-quality bases and adapter sequences from raw FASTQ reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic).
* [fastq-scan](https://github.com/rpetit3/fastq-scan) is used to estimate genome coverage of FASTQ files.
* Assembles trimmed reads into contigs using [SPAdes](https://github.com/ablab/spades).
<!--* [Kraken2](https://github.com/DerrickWood/kraken2) is used to provide detailed reports on the taxonomic composition of the trimmed raw reads. -->
* The assembled contigs from [SPAdes](https://github.com/ablab/spades) is then passed through [Prokka](https://github.com/tseemann/prokka) for annotation, [QUAST](https://quast.sourceforge.net/) for assembly statistics, [MLST](https://github.com/tseemann/mlst) for determining sequence type based on sequences of housekeeping genes, [AMRFinderPlus](https://github.com/ncbi/amr) to identify antimicrobial resistance genes, [skani](https://github.com/bluenote-1577/skani) to identify closest reference genome and [BUSCO](https://busco.ezlab.org/) for assembly completeness statistics.
* [Multiqc](https://github.com/MultiQC/MultiQC) aggregates the final outputs from [Fastqc](https://github.com/s-andrews/FastQC) , [Prokka](https://github.com/tseemann/prokka) and [QUAST](https://quast.sourceforge.net/) to produce a HTML report

The workflow generates all the output in the output prefix folder set in the config file (instructions on setup found [below](#config)). Each workflow steps gets its own individual folder as shown. **Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and some of their contents.**

```
results/2025-04-29_Project_MRSA_USA_300_raw_reads_pubQCD/
├── 2025-04-29_Project_MRSA_USA_300_raw_reads_pubQCD_Report
│   ├── data
│   │   ├── 2025-04-29_Project_MRSA_USA_300_raw_reads_pubQCD_Final_Coverage.txt
│   │   ├── 2025-04-29_Project_MRSA_USA_300_raw_reads_pubQCD_MLST_results.csv
│   │   ├── 2025-04-29_Project_MRSA_USA_300_raw_reads_pubQCD_QC_summary.csv
│   │   └── 2025-04-29_Project_MRSA_USA_300_raw_reads_pubQCD_Skani_report_final.csv
├── downsample
│   └── SRR25346107
│       ├── SRR25346107_R1_trim_paired.fastq.gz
│       └── SRR25346107_R2_trim_paired.fastq.gz
├── mlst
│   └── SRR25346107
│       └── report.tsv
├── prokka
│   └── SRR25346107
│       ├── SRR25346107.gbk
│       ├── SRR25346107.gff
│       ├── SRR25346107.tsv
│       └── SRR25346107.txt
├── quality_aftertrim
│   └── SRR25346107
│       ├── SRR25346107_Forward
│       │   ├── SRR25346107_R1_trim_paired_fastqc.html
│       │   └── SRR25346107_R1_trim_paired_fastqc.zip
├── quality_raw
│   └── SRR25346107
│       ├── SRR25346107_Forward
│       │   ├── SRR25346107_R1_fastqc.html
│       │   └── SRR25346107_R1_fastqc.zip
│       └── SRR25346107_Reverse
│           ├── SRR25346107_R2_fastqc.html
│           └── SRR25346107_R2_fastqc.zip
├── quast
│   └── SRR25346107
│       ├── quast.log
│       ├── report.html
│       ├── report.pdf
│       ├── report.tex
│       ├── report.tsv
│       ├── report.txt
│       ├── transposed_report.tex
│       ├── transposed_report.tsv
│       └── transposed_report.txt
├── raw_coverage
│   └── SRR25346107
│       └── SRR25346107_coverage.json
├── skani
│   └── SRR25346107
│       └── SRR25346107_skani_output.txt
├── spades
│   └── SRR25346107
│       ├── assembly_graph_after_simplification.gfa
│       ├── assembly_graph.fastg
│       ├── assembly_graph_with_scaffolds.gfa
│       ├── before_rr.fasta
│       ├── contigs.fasta
│       ├── spades.log
│       ├── SRR25346107_contigs_l1000.fasta
│       └── tmp
└── trimmomatic
    └── SRR25346107
        ├── SRR25346107_R1_trim_paired.fastq.gz
        ├── SRR25346107_R1_trim_paired.fastq.gz_fastqchk.txt
        ├── SRR25346107_R1_trim_unpaired.fastq.gz
        ├── SRR25346107_R2_trim_paired.fastq.gz
        └── SRR25346107_R2_trim_unpaired.fastq.gz
```


## Installation 


> If you are using Great Lakes HPC, ensure you are cloning the repository in your scratch directory. Change `your_uniqname` to your uniqname. 

```

cd /scratch/esnitkin_root/esnitkin1/your_uniqname/

```

> Clone the github directory onto your system. 

```

git clone https://github.com/Snitkin-Lab-Umich/pubQCD.git

```

> Ensure you have successfully cloned pubQCD. Type `ls` and you should see the newly created directory **_pubQCD_**. Move to the newly created directory.

```

cd pubQCD

```

> Load Bioinformatics, snakemake and singularity modules from Great Lakes modules.

```

module load Bioinformatics snakemake singularity

```
<!--
```

module load snakemake singularity

```
-->

This workflow makes use of singularity containers available through [State Public Health Bioinformatics group](https://github.com/StaPH-B/docker-builds). If you are working on Great Lakes (umich cluster)—you can load snakemake and singularity modules as shown above. However, if you are running it on your local or other computing platform, ensure you have snakemake and singularity installed.


## Setup config, cluster and SRA files

**_If you are just testing this pipeline, the config and sample files are already loaded with test data, so you do not need to make any additional changes to them. However, it is a good idea to change the prefix (name of your output folder) in the config file to give you an idea of what variables need to be modified when running your own samples on pubQCD._**

### Config
As an input, the snakemake file takes a config file where you can set the path to `sample_assembly.csv`, path to your raw sequencing reads, path to adapter fasta file etc. Instructions on how to modify `config/config.yaml` is found in `config.yaml`. 


### Cluster file

Increase/reduce the walltime depending on the number of samples you have in `config/cluster.json` to ensure the jobs are being submitted in a timely manner. 

### SRA
Download a text file of the SRA IDs of the genomes you want to run through the pipeline from NCBI/SRA/etc. and create a column called SRA_ID with the values being the sra ids.

## Quick start

### Run pubQCD on a set of samples.

> Preview the first two steps in pubQCD by performing a dryrun of the pipeline. 

```

snakemake -s workflow/download_genomes.smk --dryrun -p

```


> Submit pubQCD as a batch job. (recommended)

Change these `SBATCH` commands: `--job-name` to a more descriptive name like run_pubQCD, `--mail-user` to your email address, `--time` depending on the number of samples you have (should be more than what you specified in `cluster.json`). Feel free to make changes to the other flags if you are comfortable doing so. Once you have made the necessary changes, save the below script as `bash_script_to_run_raw_reads_pubQCD.sbat` or you can make changes directly in the slurm script in the pubQCD folder. Don't forget to submit pubQCD to Slurm! `sbatch bash_script_to_run_raw_reads_pubQCD.sbat`.

```
#!/bin/bash

#SBATCH --job-name=run_pubQCD
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@umich.edu
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=08:15:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard

# Load necessary modules
module load Bioinformatics
module load snakemake singularity

snakemake -s workflow/download_genomes.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 


# Extract prefix from the YAML config file
PREFIX=$(grep '^prefix:' config/config.yaml | awk '{print $2}')

# Define the file paths dynamically
PASS_COV_FILE="results/${PREFIX}/sample_files/samples_passed_coverage.csv"
PASS_ASSEMBLY_FILE="results/${PREFIX}/sample_files/samples_passed_assembly.csv"

echo "PASS_COV_FILE: $PASS_COV_FILE"
echo "PASS_ASSEMBLY_FILE: $PASS_ASSEMBLY_FILE"

# Run Snakemake the first time -- until coverage
snakemake -s workflow/pubQCD.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 

# samples_passed_coverage.csv should have been created in the Snakemake command above
# If not found, throw an error and exit
if [ ! -s "$PASS_COV_FILE" ] || [ "$(wc -l < "$PASS_COV_FILE")" -le 1 ]; then
    echo "Error: $PASS_COV_FILE is missing or does not have any samples. Exiting."
    exit 1
else
    echo "$PASS_COV_FILE detected. Running second part of the workflow."
fi

# Run Snakemake again to run the second part of the workflow --until assembly
snakemake -s workflow/pubQCD.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 

# samples_passed_assembly.csv should have been created in the Snakemake command above
# If not found, throw an error and exit
if [ ! -s "$PASS_ASSEMBLY_FILE" ] || [ "$(wc -l < "$PASS_ASSEMBLY_FILE")" -le 1 ]; then
    echo "Error: $PASS_ASSEMBLY_FILE is missing or does not have any samples. Exiting."
    exit 1
else
    echo "$PASS_ASSEMBLY_FILE detected. Running second part of the workflow."
fi

# Run Snakemake to finish running the rest of the pipeline
snakemake -s workflow/pubQCD.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 

# Run Snakemake for the last time to generate QC report 
snakemake -s workflow/pubQCD_report.smk -p --use-singularity --cores all

```

<!--![Alt text](./QCD_dag.svg)

### Gather Summary files and generate a report. 

>Start an interactive session in your current directory i.e. `QCD`.

```

srun --account=esnitkin1 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=5GB --cpus-per-task=1 --time=12:00:00 --pty /bin/bash

```

> Preview the steps in QCD report by performing a dryrun of the pipeline. 

```

snakemake -s QCD_report.smk --dryrun -p

```
> Run QCD report on Great lakes HPC

```

snakemake -s QCD_report.smk -p --use-singularity --cores 2

```
![Alt text](./QCD_report_dag.svg)
-->
## Dependencies

### Near Essential
* [Snakemake>=7.32.4](https://snakemake.readthedocs.io/en/stable/#)
* [Conda](https://docs.conda.io/en/latest/)

<!--All the necessary software stack required for the workflow will be installed using conda package manager.-->

### Tool stack used in workflow

* [fastq-scan](https://github.com/rpetit3/fastq-scan)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [SPades](https://github.com/ablab/spades)
<!--* [AMRFinderPlus](https://github.com/ncbi/amr)-->
* [bioawk](https://github.com/lh3/bioawk)
* [Prokka](https://github.com/tseemann/prokka)
* [mlst](https://github.com/tseemann/mlst)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
