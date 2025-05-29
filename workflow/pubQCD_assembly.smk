# Author: Ali Pirani and Dhatri Badri
configfile: "config/config_assembly.yaml"

include: "pubQCD_assembly_report.smk"

import pandas as pd
import os
import json
import glob

PREFIX = config["prefix"]

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

if not os.path.exists(f"results/"):
    os.makedirs(f"results/")


rule all:
    input:
        # spades_l1000_assembly = expand("results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta", sample=SAMPLE, prefix=PREFIX),
        prokka_gff = expand("results/{prefix}/prokka/{sample}/{sample}.gff", sample=SAMPLE, prefix=PREFIX),
        quast_report = expand("results/{prefix}/quast/{sample}/report.tsv", sample=SAMPLE, prefix=PREFIX),
        # sample_mlst_report = expand("results/{prefix}/mlst/{sample}/report.tsv", sample=SAMPLE, prefix=PREFIX),
        # skani_ref_genome_results = expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", sample=SAMPLE, prefix=PREFIX),
        # skani_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Skani_report_final.csv", prefix=PREFIX),
        # multiqc_report = expand("results/{prefix}/{prefix}_Report/multiqc/{prefix}_QC_report.html", prefix=PREFIX),
        # mlst_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_MLST_results.csv", prefix=PREFIX),
        QC_summary = expand("results/{prefix}/{prefix}_Report/data/{prefix}_QC_summary.csv", prefix=PREFIX),
 
rule bioawk:
    input:
        spades_assembly = lambda wildcards: f"{config['assembly']}/{wildcards.sample}.{config['read_extension']}"
        #spades_assembly = lambda wildcards: expand(str(config["assembly"] + "/" + f"{wildcards.sample}.{wildcards.read_extension}")),
        #spades_assembly = "results/{prefix}/spades/{sample}/{sample}.fasta"
    output:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta"
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        bioawk_params = config["bioawk"],
        prefix = "{sample}"
    # conda:
    #     "envs/bioawk.yaml"
    singularity:
        "docker://lbmc/bioawk:1.0"
    shell:
        """
        workflow/scripts/bioawk.sh {input.spades_assembly} {output.spades_l1000_assembly} {params.out_dir} {params.prefix}
        """

rule prokka:
    input:
        spades_l1000_assembly ="results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        prokka_gff = "results/{prefix}/prokka/{sample}/{sample}.gff",
    params: 
        prokka_params = config["prokka"],
        outdir = "results/{prefix}/prokka/{sample}",
        prefix = "{sample}",
    #conda:
    #    "envs/prokka.yaml"
    singularity:
        "docker://staphb/prokka:1.14.6"
    #envmodules:
    #    "Bioinformatics",
    #    "prokka"
    shell:
        "prokka -outdir {params.outdir} --strain {params.prefix} --prefix {params.prefix} {params.prokka_params} {input.spades_l1000_assembly}"

rule quast:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        quast_report = "results/{prefix}/quast/{sample}/report.tsv",
    params: 
        outdir = "results/{prefix}/quast/{sample}/",
        prefix = "{sample}",
    #conda:
    #    "envs/quast.yaml"
    singularity:
        "docker://staphb/quast:5.0.2"
    #envmodules:
    #    "Bioinformatics",
    #    "quast"
    shell:
        """
        workflow/scripts/quast.sh {input.spades_l1000_assembly} {params.outdir} 
        """
    
rule mlst:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        mlst_report = "results/{prefix}/mlst/{sample}/report.tsv",
    params: 
        outdir = "results/{prefix}/mlst/{sample}/",
        prefix = "{sample}",
    #conda:
    #    "envs/mlst.yaml"
    singularity:
        "docker://staphb/mlst:2.23.0-2024-03"
    #envmodules:
    #    "Bioinformatics",
    #    "mlst"
    shell:
        "mlst {input.spades_l1000_assembly} > {output.mlst_report}"

rule busco:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        busco_out = "results/{prefix}/{sample}/busco/busco.txt",
    params: 
        outdir = "results/{prefix}/busco/{sample}/",
        prefix = "{sample}",
        threads = config["ncores"],
    #conda:
    #    "envs/busco.yaml"
    singularity:
        "docker://staphb/busco:5.7.1-prok-bacteria_odb10_2024-01-08"
    #envmodules:
    #    "Bioinformatics",
    #    "busco"
    shell:
        "busco -f -i {input.spades_l1000_assembly} -m genome -l bacteria_odb10 -o {params.outdir}"

rule skani:
    input:
        spades_assembly = lambda wildcards: f"{config['assembly']}/{wildcards.sample}.{config['read_extension']}"
        #spades_assembly = lambda wildcards: expand(str(config["assembly"] + "/" + f"{wildcards.sample}.fasta")),
    output:
        skani_output = "results/{prefix}/skani/{sample}/{sample}_skani_output.txt"
    params:
        skani_ani_db = config["skani_db"],
        threads = 6
    #conda:
    #    "envs/skani.yaml"
    singularity:
        "docker://staphb/skani:0.2.1"
    shell:
        "skani search {input.spades_assembly} -d {params.skani_ani_db} -o {output.skani_output} -t {params.threads}"



"""
END OF PIPELINE
"""