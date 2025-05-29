# Author: Dhatri Badri
configfile: "config/config.yaml"

import pandas as pd
import os
import json

samples_df = pd.read_csv(config["download_genomes"])
SRA_IDS = list(samples_df['SRA_ID'])

rule all:
    input:
        # expand(os.path.join(config["downloaded_genomes_outdir"], "{sra_id}", "{sra_id}.sra"), sra_id=SRA_IDS),
        expand(os.path.join(config["downloaded_genomes_outdir"], "{sra_id}_R1.fastq.gz"), sra_id=SRA_IDS),

rule download_sra:
    output:
        temp(os.path.join(config["downloaded_genomes_outdir"], "{sra_id}", "{sra_id}.sra"))
    params:
        outdir = config["downloaded_genomes_outdir"],
        sra_id = lambda wildcards: wildcards.sra_id
    singularity:
        "docker://staphb/sratoolkit:3.2.1"
    shell:
        """
        cd {params.outdir}
        prefetch -O {params.outdir} {params.sra_id}
        """

rule fastq_dump:
    input:
        sra=os.path.join(config["downloaded_genomes_outdir"], "{sra_id}","{sra_id}.sra")
    output:
        R1=os.path.join(config["downloaded_genomes_outdir"], "{sra_id}_R1.fastq.gz")
    params:
        outdir=config["downloaded_genomes_outdir"],
        sra_id = lambda wildcards: wildcards.sra_id
    singularity:
        "docker://staphb/sratoolkit:3.2.1"
    shell:
        """
        cd {params.outdir}
        fastq-dump --outdir {params.outdir} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {input.sra}
        mv {params.sra_id}_pass_1.fastq.gz {params.sra_id}_R1.fastq.gz
        mv {params.sra_id}_pass_2.fastq.gz {params.sra_id}_R2.fastq.gz
        """

"""
END OF PIPELINE
"""