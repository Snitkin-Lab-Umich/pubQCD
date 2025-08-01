# Author: Ali Pirani and Dhatri Badri  
#configfile: "config/config_assembly.yaml"

import pandas as pd
import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re

PREFIX = config["prefix"]

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

# if not os.path.exists("results/" + PREFIX):
    # os.system("mkdir %s" % "results/" + PREFIX)

# Organize reports directory
prefix = PREFIX
outdir = "results/%s" % prefix
report_dir = outdir + "/%s_Report" % prefix
report_script_dir = report_dir + "/scripts"
report_data_dir = report_dir + "/data"
report_multiqc_dir = report_dir + "/multiqc"
report_fig_dir = report_dir + "/fig"

isExist = os.path.exists(report_dir)
if not isExist:
    os.makedirs(report_dir)

isExist = os.path.exists(report_script_dir)
if not isExist:
    os.makedirs(report_script_dir)

isExist = os.path.exists(report_data_dir)
if not isExist:
    os.makedirs(report_data_dir)

isExist = os.path.exists(report_multiqc_dir)
if not isExist:
    os.makedirs(report_multiqc_dir)

isExist = os.path.exists(report_fig_dir)
if not isExist:
    os.makedirs(report_fig_dir)

# def coverage_report(prefix, outdir):
#     prefix = prefix.pop()
#     report_dir = str(outdir.pop()) + "/%s_Report" % prefix
#     # Generate Coverage report 
#     final_coverage_file = "%s/data/%s_Final_Coverage.txt" % (report_dir, prefix)
#     f3=open(final_coverage_file, 'w+')
#     header = "Sample,Total_reads,Total_bp,MeanReadLength,Coverage\n"
#     f3.write(header)

#     for sampl in SAMPLE_PASS_ASSEMBLY:
#         coverage_json = "results/%s/raw_coverage/%s/%s_coverage.json" % (prefix, sampl, sampl)
#         f = open(coverage_json)
#         data = json.load(f)
#         # data = json.loads(coverage_json)
#         f3.write("%s,%s,%s,%s,%s\n" % (sampl, data['qc_stats']['read_total'], data['qc_stats']['total_bp'], data['qc_stats']['read_mean'], data['qc_stats']['coverage']))
#     f3.close()   

#     Coverage = pd.read_csv(final_coverage_file, sep=',', header=0)
#     Coverage = Coverage.replace(['_R1.fastq.gz'], '', regex=True)

    #print ("Number of Samples in Coverage Report - %s" % len(Coverage))

    #Coverage_NEG_CNTL = Coverage[Coverage.Sample.str.match('(.*NEG*)')]

    #print ("Number of Negative Control samples %s" % len(Coverage_NEG_CNTL))

    #print ("Number of Negative Control samples with > 100X coverage %s" % len(Coverage_NEG_CNTL[Coverage_NEG_CNTL['Coverage'] > 100]))

    #Coverage_dist = Coverage.sort_values(by='Coverage',ascending=False).plot(x='Sample_name', y='Coverage', kind="barh", title="Estimated Genome Coverage", figsize=(20, 20), fontsize=40).get_figure()

    #Coverage_dist.savefig('%s/%s_Coverage_distribution.pdf' % (report_dir, prefix))

#def kraken_report(prefix, outdir):
#    prefix = prefix.pop()
#    outdir = outdir.pop()
    
    # Organize reports directory
#    report_dir = str(outdir) + "/%s_Report" % prefix
#    report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix

#    kraken_dir = str(outdir) + "/*/kraken"

#    kraken_summary_script = open("%s/kraken_summary.sh" % report_script_dir, 'w+')
#    kraken_summary_script.write("echo \"Sample,Percentage of reads for Species,# of reads for Species, Species\" > %s/data/%s_Kraken_report_final.csv\n" % (report_dir, prefix))
#    kraken_summary_script.write("for i in results/%s/*/kraken/*_kraken_report.tsv; do grep -w 'S' $i | sort -k1n | tail -n1; done > /tmp/Kraken_report_temp.txt\n" % prefix)
#    kraken_summary_script.write("ls -d results/%s/*/kraken/*_kraken_report.tsv | awk -F'/' '{print $NF}' | sed 's/_kraken_report.tsv//g' > %s/data/samplenames.txt\n" % (prefix, report_dir))
#    kraken_summary_script.write("paste %s/data/samplenames.txt /tmp/Kraken_report_temp.txt > /tmp/Kraken_report_combined.txt\n" % (report_dir))
#    kraken_summary_script.write("awk -F'\\t' 'BEGIN{OFS=\",\"};{print $1, $2, $3, $7}' /tmp/Kraken_report_combined.txt >> %s/data/%s_Kraken_report_final.csv\n" % (report_dir, prefix))
#    kraken_summary_script.write("sed -i 's/\s//g' %s/data/%s_Kraken_report_final.csv\n" % (report_dir, prefix))
#    kraken_summary_script.close()

#    os.system("bash %s/kraken_summary.sh" % report_script_dir)

def skani_report(outdir, prefix):
    prefix = prefix.pop()
    outdir = "results/%s" % prefix
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_data_dir = report_dir + "/data"
    result_df = pd.DataFrame(columns=['Sample', 'ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name', 'Species'])  # Add 'Species' column

    skani_dir = os.path.join(outdir, 'skani')  # Navigate to skani directory

    for sample_name in os.listdir(skani_dir):  # Iterate over samples in the results/prefix/skani directory
        sample_dir = os.path.join(skani_dir, sample_name)

        if os.path.isdir(sample_dir):  # Check if it's a directory
            skani_file_path = os.path.join(sample_dir, f'{sample_name}_skani_output.txt')  # Look for the skani output file

            if os.path.exists(skani_file_path):  # Check if the skani file exists
                skani_file = pd.read_csv(skani_file_path, sep='\t| ,', skipinitialspace=True, header=0)  # Read the skani file
                first_row_df = skani_file[['ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name']].iloc[:1]  # Extract the first row

                if first_row_df.empty:  # Check if the first row is empty
                    first_row_df = pd.DataFrame({
                        'Sample': [sample_name],  # Add sample name
                        'ANI': ["NA"],
                        'Align_fraction_ref': ["NA"],
                        'Align_fraction_query': ["NA"],
                        'Ref_name': ["NA"],
                        'Species': ["NA"]  # Add NAs for Species
                    })
                else:
                    first_row_df.loc[:, 'Sample'] = sample_name  # Add sample name
                    # Extract species using regex from Ref_name
                    first_row_df.loc[:, 'Species'] = first_row_df['Ref_name'].apply(
                        lambda x: re.search(r"[A-Za-z]+\s[A-Za-z]+", x).group(0) if pd.notnull(x) and re.search(r"[A-Za-z]+\s[A-Za-z]+", x) else "NAs"
                    )

                first_row_df = first_row_df[['Sample', 'ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name', 'Species']]  # Reorder columns
                result_df = pd.concat([result_df, first_row_df], ignore_index=True)  # Concatenate to the result dataframe

    result_file_path = os.path.join(report_data_dir, f'{prefix}_Skani_report_final.csv')  # Save final result to CSV
    result_df.to_csv(result_file_path, index=False)


def summary(prefix, outdir, skani_genome_size):
    prefix = prefix.pop()
    outdir = outdir.pop()
    
    # Organize reports directory
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix
    
    
    # Coverage = pd.read_csv("results/%s/%s_Report/data/%s_Final_Coverage.txt" % (prefix, prefix, prefix), sep=',', header=0)
    # Coverage.rename(columns = {'Sample_name':'Sample'}, inplace = True)

    #kraken = pd.read_csv("results/%s/%s_Report/data/%s_Kraken_report_final.csv" % (prefix, prefix, prefix), sep=',', header=0)
    
    mlst = pd.read_csv("results/%s/%s_Report/data/%s_MLST_results.csv" % (prefix, prefix, prefix), sep='\t', header=0)
    #mlst = mlst.replace(['_contigs_l1000.fasta'], '', regex=True)
    #mlst = mlst.replace(['results/.*/spades/'], '', regex=True)
    #mlst = mlst.replace(['%s' % prefix], '', regex=True)
    mlst['Sample'] = mlst['Sample'].replace(r'.*/spades/(.*?)/.*', r'\1', regex=True)

    # multiqc_fastqc_summary = pd.read_csv("results/%s/%s_Report/multiqc/%s_QC_report_data/multiqc_fastqc.txt" % (prefix, prefix, prefix), sep='\t', header=0)
    # patternDel = "_R2"
    # filter = multiqc_fastqc_summary['Sample'].str.contains(patternDel)
    # multiqc_fastqc_summary = multiqc_fastqc_summary[~filter]
    # aftertrim_filter = multiqc_fastqc_summary['Sample'].str.contains("_R1_trim_paired")
    # raw_multiqc_fastqc_summary = multiqc_fastqc_summary[~aftertrim_filter]
    # raw_multiqc_fastqc_summary = raw_multiqc_fastqc_summary.replace(['_R1'], '', regex=True)
    
    # aftertrim_multiqc_fastqc_summary = multiqc_fastqc_summary[aftertrim_filter]
    # aftertrim_multiqc_fastqc_summary = aftertrim_multiqc_fastqc_summary.replace(['_R1_trim_paired'], '', regex=True)
    # aftertrim_multiqc_fastqc_summary = aftertrim_multiqc_fastqc_summary.add_prefix('After_trim_')
    # aftertrim_multiqc_fastqc_summary.rename(columns = {'After_trim_Sample':'Sample'}, inplace = True)

    multiqc_general_stats_summary = pd.read_csv("results/%s/%s_Report/multiqc/%s_QC_report_data/multiqc_general_stats.txt" % (prefix, prefix, prefix), sep='\t', header=0)
    quast_filter = multiqc_general_stats_summary['Sample'].str.contains("_contigs_l1000")
    multiqc_quast = multiqc_general_stats_summary[quast_filter]
    multiqc_quast = multiqc_quast.replace(['_contigs_l1000'], '', regex=True)
    
    if 'QUAST_mqc-generalstats-quast-N50' in multiqc_quast.columns and 'QUAST_mqc-generalstats-quast-Total_length' in multiqc_quast.columns:
        multiqc_quast = multiqc_quast[["Sample", "QUAST_mqc-generalstats-quast-N50", "QUAST_mqc-generalstats-quast-Total_length"]]
        multiqc_quast = multiqc_quast.rename(columns={"QUAST_mqc-generalstats-quast-N50": "N50", "QUAST_mqc-generalstats-quast-Total_length": "Total length"})
    elif 'N50' in multiqc_quast.columns and 'Total length' in multiqc_quast.columns:
        multiqc_quast = multiqc_quast[["Sample", "N50", "Total length"]]
    #multiqc_quast = multiqc_quast[["Sample", "N50", "Total length"]]

    contig_distribution = pd.read_csv("results/%s/%s_Report/multiqc/%s_QC_report_data/mqc_quast_num_contigs_1.txt" % (prefix, prefix, prefix), sep='\t', header=0)
    contig_distribution = contig_distribution.replace(['_contigs_l1000'], '', regex=True)
    contig_distribution['Total # of contigs'] = contig_distribution.sum(axis=1, numeric_only=True)
    contig_distribution = contig_distribution[['Sample', 'Total # of contigs']]
    
    #read final skani output file
    skani_summary = pd.read_csv("results/%s/%s_Report/data/%s_Skani_report_final.csv" % (prefix, prefix, prefix), sep=',', skipinitialspace=True, header=0, engine='python')

    # QC_summary_temp1 = pd.merge(Coverage, mlst, on=["Sample", "Sample"],  how='left')
    # QC_summary_temp2 = QC_summary_temp1
    # QC_summary_temp3 = pd.merge(QC_summary_temp2, raw_multiqc_fastqc_summary, on=["Sample", "Sample"], how='left')
    # QC_summary_temp4 = pd.merge(QC_summary_temp3, aftertrim_multiqc_fastqc_summary, on=["Sample", "Sample"], how='left')
    QC_summary_temp1 = pd.merge(mlst, multiqc_quast, on=["Sample", "Sample"], how='left')
    QC_summary_temp2 = pd.merge(QC_summary_temp1, contig_distribution, on=["Sample", "Sample"], how='left')
    
    #QC_summary_temp7 = QC_summary_temp6[["Sample" , "Total_reads" , "Total_bp" , "MeanReadLength" , "Coverage" , "Scheme" , "ST" , "PercentageofreadsforSpecies" , "#ofreadsforSpecies" , "Species" , "After_trim_per_base_sequence_content" , "After_trim_overrepresented_sequences" , "After_trim_%GC" , "After_trim_Total Bases" , "After_trim_Total Sequences" , "After_trim_median_sequence_length" , "After_trim_avg_sequence_length" , "After_trim_total_deduplicated_percentage" , "After_trim_Sequence length" , "After_trim_adapter_content" , "N50" , "Total length" , "Total # of contigs"]].copy() #.copy() to deal with SettingWithCopyWarning error
    QC_summary_temp3 = QC_summary_temp2[["Sample" , "Scheme" , "ST" , "N50" , "Total length" , "Total # of contigs"]].copy() #.copy() to deal with SettingWithCopyWarning error

    skani_genome_size = list(skani_genome_size)[0]
    
    # Read in skani species gneome size table
    skani_genome_table = pd.read_csv(skani_genome_size)

    # Check the total length against the assembly length with ±15% rule
    def check_assembly_length(total_length, assembly_length):
        if pd.isnull(total_length):
            return 'FAIL'
        if pd.isnull(assembly_length):
            if config["genome_size"] <= total_length <= config["assembly_length"]:
                return 'PASS'
            else:
                return 'FAIL'
        lower_bound = assembly_length * 0.85
        upper_bound = assembly_length * 1.15
        if lower_bound <= total_length <= upper_bound:
            return 'PASS'
        else:
            return 'FAIL'
    
    QC_summary_temp4 = pd.merge(QC_summary_temp3, skani_summary, on=["Sample", "Sample"], how='left') # Merge skani df into the existing dataframe

    # Merge QC summary with skani_genome_size on species
    QC_summary_temp4 = QC_summary_temp4.merge(
        skani_genome_table,
        left_on='Species',  
        right_on='Species', 
        how='left'
    )

    # Apply the length check function
    QC_summary_temp4['Length Check'] = QC_summary_temp4.apply(
        lambda row: check_assembly_length(row['Total length'], row['Assembly_Length']),
        axis=1
    )

    # Updated QC check based on the new Length Check condition
    QC_check_condition = [
        (QC_summary_temp4['Total # of contigs'] > config["max_contigs"]),
        (QC_summary_temp4['Total # of contigs'] < config["min_contigs"]),
        #(QC_summary_temp4['Coverage'] < config["coverage"]),
        (QC_summary_temp4['Length Check'] == 'FAIL'),
        (QC_summary_temp4['Total # of contigs'].isnull()),
        (pd.isnull(QC_summary_temp4['Total length'])),
    ]

    status = ['FAIL', 'FAIL', 'FAIL', 'FAIL', "Run FAIL"]

    QC_summary_temp4['QC Check'] = np.select(QC_check_condition, status, default='PASS')
    
     # Remove the 'Length Check' column
    QC_summary_temp4 = QC_summary_temp4.drop(columns=['Length Check', 'Assembly_Length'])

    # First, get the current list of columns
    columns = list(QC_summary_temp4.columns)
    
    # Insert QC Check between Total # of contigs and ANIs
    # Need to know the index positions of these columns to make the rearrangement correctly
    contigs_index = columns.index('Total # of contigs')
    ani_index = columns.index('ANI')
    qc_check_index = columns.index('QC Check')

    # Create the new column order
    # Put all columns before Total # of contigs then Total # of contigs, QC Check, ANI and the rest
    new_columns = columns[:contigs_index + 1] + ['QC Check'] + columns[ani_index:qc_check_index]

    # Rearrange the columns
    QC_summary_temp4 = QC_summary_temp4[new_columns]

    QC_summary_temp4.to_csv('results/%s/%s_Report/data/%s_QC_summary.csv' % (prefix, prefix, prefix), index=False)

def plot(prefix, outdir):
    prefix = prefix.pop()
    outdir = outdir.pop()
    
    # Organize reports directory
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix
    
    QC_summary = pd.read_csv('results/%s/%s_Report/data/%s_QC_summary.csv' % (prefix, prefix, prefix), sep=',', header=0)    

    Coverage = pd.read_csv("results/%s/%s_Report/data/%s_Final_Coverage.txt" % (prefix, prefix, prefix), sep=',', header=0)    
    Coverage_dist = QC_summary.sort_values(by='Coverage',ascending=False).plot(x='Sample', y='Coverage', kind="barh", title="Estimated Genome Coverage", figsize=(20, 20), fontsize=40).get_figure()
    Coverage_dist.savefig('%s/fig/%s_Coverage_distribution.png' % (report_dir, prefix), dpi=600)


    ax1 = QC_summary.plot.scatter(x = 'After_trim_total_deduplicated_percentage', y = 'After_trim_Total Sequences', c = 'DarkBlue')
    fig = ax1.get_figure()
    fig.savefig('%s/fig/%s_raw_dedup_vs_totalsequence.png' % (report_dir, prefix), dpi=600)

    ax1 = QC_summary.plot.scatter(x = 'After_trim_total_deduplicated_percentage', y = 'After_trim_Total Sequences', c = 'DarkBlue')
    fig = ax1.get_figure()
    fig.savefig('%s/fig/%s_aftertrim_dedup_vs_totalsequence.png' % (report_dir, prefix), dpi=600)
    ax1.cla()

    #ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['After_trim_%GC'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=600)
    #plt.savefig('%s/fig/%s_Assembly_contig_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=200)
    #ax.cla()

    #ax = sns.scatterplot(x=QC_summary['Total length'], y=QC_summary['After_trim_%GC'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=600)
    #plt.savefig('%s/fig/%s_Assembly_length_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=200)
    #ax.cla()

    #ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['N50'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=600)
    #plt.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=200)
    #ax.cla()

    #ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['Coverage'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=600)
    #plt.savefig('%s/fig/%s_Assembly_contig_vs_Coverage.png' % (report_dir, prefix), dpi=200)
    #ax.cla()

    #ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['Total length'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=600)
    #plt.savefig('%s/fig/%s_Assembly_contig_vs_length.png' % (report_dir, prefix), dpi=200)
    #ax.cla()

# rule all:
#     input:
#         # coverage_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Final_Coverage.txt", prefix=PREFIX),
#         #kraken_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Kraken_report_final.csv", prefix=PREFIX),
#         skani_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Skani_report_final.csv", prefix=PREFIX),
#         multiqc_report = expand("results/{prefix}/{prefix}_Report/multiqc/{prefix}_QC_report.html", prefix=PREFIX),
#         mlst_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_MLST_results.csv", prefix=PREFIX),
#         QC_summary = expand("results/{prefix}/{prefix}_Report/data/{prefix}_QC_summary.csv", prefix=PREFIX),
#         #QC_plot = expand("results/{prefix}/{prefix}_Report/fig/{prefix}_Coverage_distribution.png", prefix=PREFIX)

# rule coverage_report:
#     input:
#         outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
#         coverage_out = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", prefix=PREFIX, sample=SAMPLE_PASS_ASSEMBLY)
#     output:
#         coverage = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Final_Coverage.txt",
#     params:
#         prefix = "{prefix}",
#     run:
#         coverage_report({params.prefix}, {input.outdir})

# rule amr_report:
#     input:
#         outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
#     output:
#         amr_summary = f"results/{{prefix}}/report/{{prefix}}_AMR_minimal_report.csv",
#     params:
#         prefix = "{prefix}",
#         phandango = "--no_tree"
#     conda:
#         "envs/ariba.yaml"
#     #singularity:
#     #    "docker://staphb/ariba:2.14.7"
#     shell:
#         "ariba summary --preset minimal {params.phandango} {input.outdir}/report/{params.prefix}_AMR_minimal_report {input.outdir}/*/ariba_card/report.tsv && ariba summary --preset all {params.phandango} {input.outdir}/report/{params.prefix}_AMR_all_report {input.outdir}/*/ariba_card/report.tsv"

#rule kraken_report:
#    input:
#        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
#    output:
#        kraken_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Kraken_report_final.csv",
#    params:
#        prefix = "{prefix}",
#    run:
#        kraken_report({params.prefix}, {input.outdir})

rule skani_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        skani_out = expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", prefix=PREFIX, sample=SAMPLE)
    output:
        skani_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Skani_report_final.csv",
    params:
        prefix = "{prefix}",
    run:
        skani_report({input.outdir}, {params.prefix})

rule multiqc:
    input:
        inputdir = lambda wildcards: expand(f"results/{wildcards.prefix}"),
        # coverage = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Final_Coverage.txt"),
        #kraken = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Kraken_report_final.csv"),
        mlst = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_MLST_results.csv"),
    output:
        multiqc_fastqc_report = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report.html",
        # multiqc_fastqc = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_fastqc.txt",
        # multiqc_general_stats = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_general_stats.txt",
        #fastqc_report = f"results//{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_fastqc.txt"
    params:
        outdir = "results/{prefix}/{prefix}_Report",
        prefix = "{prefix}",
    #conda:
    #    "envs/multiqc.yaml"
    singularity:
        "docker://staphb/multiqc:1.19"
    shell:
        "multiqc -f -s --export --outdir {params.outdir}/multiqc -n {params.prefix}_QC_report -i {params.prefix}_QC_report {input.inputdir}/prokka/* {input.inputdir}/quast/*"

rule mlst_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        mlst_out = expand("results/{prefix}/mlst/{sample}/report.tsv", prefix=PREFIX, sample=SAMPLE)
    output:
        mlst_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_MLST_results.csv",
    params:
        prefix = "{prefix}",
    shell:
        "echo \"Sample\tScheme\tST\" > {output.mlst_report} && cut -f1-3 {input.outdir}/mlst/*/report.tsv >> {output.mlst_report}"

rule Summary:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        multiqc_fastqc_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/multiqc/{wildcards.prefix}_QC_report.html"),
        # multiqc_fastqc = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/multiqc/{wildcards.prefix}_QC_report_data/multiqc_fastqc.txt"),
        # multiqc_general_stats = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/multiqc/{wildcards.prefix}_QC_report_data/multiqc_general_stats.txt"),
        # coverage = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Final_Coverage.txt"),
        #kraken = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Kraken_report_final.csv"),
        mlst = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_MLST_results.csv"),
        skani_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Skani_report_final.csv"),
    output:
        QC_summary_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_QC_summary.csv",
    params:
        prefix = "{prefix}",
        skani_genome_size_table = config["skani_genome_size"] # In the config folder
    run:
        summary({params.prefix}, {input.outdir}, {params.skani_genome_size_table})

rule plot:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        QC_summary_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_QC_summary.csv"),
    output:
        QC_summary_report = f"results/{{prefix}}/{{prefix}}_Report/fig/{{prefix}}_Coverage_distribution.png",
    params:
        prefix = "{prefix}",
    run:
        plot({params.prefix}, {input.outdir})
