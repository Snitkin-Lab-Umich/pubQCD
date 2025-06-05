# pubQCD - Public datasets Quality Control and Contamination Detection workflow

pubQCD is a quality control pipeline for datasets downloaded from public database implemented in Snakemake that takes either raw fastqs or assemblies as input.

If you want to process raw reads, click [here](README.md), otherwise, continue reading the instructions on how to process assemblies. **The pipeline assumes you have downloded the assemblies beforehand and are pointing to the directory where these genomes are.**

### Summary

As part of the SOP in the [Snitkin lab](https://thesnitkinlab.com/index.php), this pipeline should be run on raw sequencing data (separate pipeline for assemblies coming soon!) as soon as the data is downloaded from a public database. 

In short, it performs the following steps:

* Assemblies are filtered for contigs >1kb using [bioawk](https://github.com/lh3/bioawk).
* The contigs from [bioawk](https://github.com/lh3/bioawk) are then passed through [Prokka](https://github.com/tseemann/prokka) for annotation, [QUAST](https://quast.sourceforge.net/) for assembly statistics, [MLST](https://github.com/tseemann/mlst) for determining sequence type based on sequences of housekeeping genes, [skani](https://github.com/bluenote-1577/skani) to identify closest reference genome and [BUSCO](https://busco.ezlab.org/) for assembly completeness statistics.
* [Multiqc](https://github.com/MultiQC/MultiQC) aggregates the final outputs from [Prokka](https://github.com/tseemann/prokka) and [QUAST](https://quast.sourceforge.net/) to produce a HTML report
* An additional QC report is generated as the final step in the pipeline that will indicate which samples have passed and failed the QC pipeline. 

The workflow generates all the output in the output prefix folder set in the config file (instructions on setup found [below](#config)). Each workflow steps gets its own individual folder as shown. **Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and some of their contents.**

```
results/2025-05-29_Project_MRSA_USA_300_assembly_pubQCD/
├── 2025-05-29_Project_MRSA_USA_300_assembly_pubQCD_Report
│   ├── data
│   │   ├── 2025-05-29_Project_MRSA_USA_300_assembly_pubQCD_MLST_results.csv
│   │   ├── 2025-05-29_Project_MRSA_USA_300_assembly_pubQCD_QC_summary.csv
│   │   └── 2025-05-29_Project_MRSA_USA_300_assembly_pubQCD_Skani_report_final.csv
├── mlst
│   └── GCA_000245595.2_ASM24559v2_genomic
│       └── report.tsv
├── prokka
│   └── GCA_000245595.2_ASM24559v2_genomic

│       ├── GCA_000245595.2_ASM24559v2_genomic.gbk
│       ├── GCA_000245595.2_ASM24559v2_genomic.gff
│       └── GCA_000245595.2_ASM24559v2_genomic.txt
├── quast
│   └── GCA_000245595.2_ASM24559v2_genomic
│       ├── quast.log
│       ├── report.html
│       ├── report.pdf
│       ├── report.tex
│       ├── report.tsv
│       ├── report.txt
│       ├── transposed_report.tex
│       ├── transposed_report.tsv
│       └── transposed_report.txt
├── skani
│   └── GCA_000245595.2_ASM24559v2_genomic
│       └── GCA_000245595.2_ASM24559v2_genomic_skani_output.txt
└── spades
    └── GCA_000245595.2_ASM24559v2_genomic
        ├── GCA_000245595.2_ASM24559v2_genomic_contigs_l1000.fasta
        └── spades_assembly_header_info.txt
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
> Ensure scripts to use in pipeline are executable
```

chmod +x workflow/scripts/bioawk.sh
chmod +x workflow/scripts/quast.sh

```


This workflow makes use of singularity containers available through [State Public Health Bioinformatics group](https://github.com/StaPH-B/docker-builds). If you are working on Great Lakes (umich cluster)—you can load snakemake and singularity modules as shown above. However, if you are running it on your local or other computing platform, ensure you have snakemake and singularity installed.


## Setup config, samples and cluster files

**_If you are just testing this pipeline, the config and sample files are already loaded with test data, so you do not need to make any additional changes to them. However, it is a good idea to change the prefix (name of your output folder) in the config file to give you an idea of what variables need to be modified when running your own samples on pubQCD._**

### Config
As an input, the snakemake file takes a config file where you can set the path to `sample_assembly.csv`, path to your raw sequencing reads, path to adapter fasta file etc. Instructions on how to modify `config/config_assembly.yaml` is found in `config_assembly.yaml`. 

### Samples
Add samples to `config/sample_assembly.csv` following the explanation provided below. `sample_assembly.csv` should be a comma seperated file consisting of two columns—`sample_id` and `illumina_r1`.

* `sample_id` is the prefix that should be extracted from your assemblies. For example, in  your assembly directory, if you have a file called `GCA_000177995.1_ASM17799v1_genomic.fna`, your sample_id would be `GCA_000177995.1_ASM17799v1_genomic`.

You can create sample_assembly.csv file using the following for loop. Replace *path_to_your_raw_reads* below with the actual path to your raw sequencing reads and *suffix* with the suffix of your assemblies either `.fna` or `.fasta`

```

echo "sample_id" > config/sample_assembly.csv

for read1 in path_to_your_raw_reads/*suffix; do
    sample_id=$(basename $read1 | sed 's/suffix//g')
    echo $sample_id
done >> config/sample_assembly.csv

```

### Cluster file

## Quick start


### Run pubQCD on a set of samples.

> Preview the first two steps in pubQCD by performing a dryrun of the pipeline. 

```

snakemake -s workflow/pubQCD_assembly.smk --dryrun -p

```

>Run pubQCD on terminal directly. 

```

snakemake -s workflow/pubQCD_assembly.smk -p --use-conda --use-singularity --use-envmodules -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=assembly_slurm_out/slurm-%j.out" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config_assembly.yaml --latency-wait 1000 --nolock 

```
> Submit pubQCD as a batch job. (recommended)

Change these `SBATCH` commands: `--job-name` to a more descriptive name like run_pubQCD, `--mail-user` to your email address, `--time` depending on the number of samples you have (should be more than what you specified in `cluster.json`). Feel free to make changes to the other flags if you are comfortable doing so. Once you have made the necessary changes, save the below script as `bash_script_to_run_assembly_pubQCD.sbat` or you can make changes directly in the slurm script in the pubQCD folder. Don't forget to submit pubQCD to Slurm! `sbatch bash_script_to_run_assembly_pubQCD.sbat`.

```
#!/bin/bash

#SBATCH --job-name=run_pubQCD
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@umich.edu
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=08:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard

# Load necessary modules
module load Bioinformatics
module load snakemake singularity

# Run main pipeline
snakemake -s workflow/pubQCD_assembly.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=assembly_slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config_assembly.yaml --latency-wait 1000 --nolock 


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
![Alt text](images/pubQCD_assembly_dag.svg)
-->
## Dependencies

### Near Essential
* [Snakemake>=7.32.4](https://snakemake.readthedocs.io/en/stable/#)
* [Conda](https://docs.conda.io/en/latest/)

<!--All the necessary software stack required for the workflow will be installed using conda package manager.-->

### Tool stack used in workflow

* [BUSCO](https://busco.ezlab.org/)
* [bioawk](https://github.com/lh3/bioawk)
* [Prokka](https://github.com/tseemann/prokka)
* [mlst](https://github.com/tseemann/mlst)
* [MultiQC](https://multiqc.info/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
