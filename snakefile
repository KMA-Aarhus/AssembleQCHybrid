# snakemake --profile configs/slurm


print("/*")
__author__ = "Tine Ebsen" # Please add your name here if you make changes.
__version__ = "0.3"

from pathlib import Path
import os
import pandas as pd
import sys





print("        /\\                          | |   | |     / __ \\ / ____| KMA,AUH")
print("       /  \\   ___ ___  ___ _ __ ___ | |__ | | ___| |  | | |      ")
print("      / /\\ \\ / __/ __|/ _ \\ '_ ` _ \\| '_ \\| |/ _ \\ |  | | |     ")
print("     / ____ \\__ \\__ \\  __/ | | | | | |_) | |  __/ |__| | |____  ")
print("    /_/    \\_\\___/___/\\___|_| |_| |_|_.__/|_|\\___|\\___\\_\\_____| ")
print(" ")


# Field variables, which might make sense to set from a config file.
#clean_uploads_dir = "../BACKUP/nanopore_sarscov2/pappenheim_clean"
#clean_uploads_dir = "../BACKUP/nanopore_sarscov2/pappenheim_clean/testdir"

out_base = config["out_base"]
sample_reads = config["sample_reads"]
samplename = config["samplename"]
R1 = config["R1"]
R2 = config["R2"]
ONT = config["ONT"]

if config["samplename"] == "NA":
    raise Exception("No samplename file was given. Please specify a samplename by appending --config samplename=\"name\" to the command line call.")
if config["R1"] == "NA":
    raise Exception("No R1 reads were given. Please specify R1 reads by appending --config rundir=\"path/to/R1 reads/\" to the command line call.")
if config["R2"] == "NA":
    raise Exception("No R2 reads were given. Please specify R2 reads by appending --config rundir=\"path/to/R2 reads/\" to the command line call.")
if config["ONT"] == "NA":
    raise Exception("No ONT reads were given. Please specify ONT reads by appending --config rundir=\"path/to/ONT reads/\" to the command line call.")


kraken2_db = config["kraken2_db"]
plasmidfinder_db = config["plasmidfinder_db"]

#df = pd.DataFrame(files, columns = ["filename"])
print([samplename, R1, R2, ONT])
df = pd.DataFrame([[samplename, R1, R2, ONT]], columns = ["sample_id", "R1", "R2", "ONT"])
print(df)

rule all:
    input:
        expand(["{out_base}/{sample_id}/trimmed/{sample_id}_val_1.fq.gz", \
                "{out_base}/{sample_id}/trimmed/{sample_id}_val_2.fq.gz", \
                "{out_base}/{sample_id}/trimmed/{sample_id}_ONT_trimmed.fq.gz", \
                "{out_base}/{sample_id}/qc/{sample_id}_nanostat", \
                "{out_base}/{sample_id}/{sample_id}_paired_kraken2_reads_report.txt", \
                "{out_base}/{sample_id}/{sample_id}_ONT_kraken2_reads_report.txt", \
                "{out_base}/{sample_id}/assembly/{sample_id}_report.txt", \
                "{out_base}/{sample_id}/{sample_id}_consensus.fasta", \
                "{out_base}/{sample_id}/{sample_id}_consensus.gff", \
                "{out_base}/{sample_id}/{sample_id}_consensus.gbk", \
                "{out_base}/{sample_id}/plasmidfinder/{sample_id}_data.json", \
                "{out_base}/multiqc_report.html" \
                ], \
                
               out_base = out_base, sample_id = df["sample_id"])





####################
# Setup for data analysis #
####################

rule nanostat:
    input:
        ONT = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["ONT"].values[0]
    output: 
        "{out_base}/{sample_id}/qc/{sample_id}_nanostat"
    conda: "configs/conda.yaml"
    threads: 1
    shell: """

    NanoStat --fastq {input.ONT} -o {out_base}/{wildcards.sample_id}/qc/ -n {wildcards.sample_id}_nanostat

    """



    # Trim adapters
rule trim_adapt_PE:
    input: 
        R1 = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["R1"].values[0],
        R2 = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["R2"].values[0]
    output: 
        R1 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_1.fq.gz",
        R2 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_2.fq.gz"
    conda: "configs/conda.yaml"
    threads: 4
    shell: """
        mkdir -p {out_base}/{wildcards.sample_id}/trimmed

        trim_galore --paired --gzip --cores 4 --basename {wildcards.sample_id} --fastqc -o {out_base}/{wildcards.sample_id}/trimmed  {input.R1} {input.R2} --length 100 --quality 25
            
        """

rule trim_adapt_ONT:
    input: 
        ONT = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["ONT"].values[0]
    output: 
        ONT = "{out_base}/{sample_id}/trimmed/{sample_id}_ONT_trimmed.fq.gz"
    conda: "configs/conda.yaml"
    threads: 4
    shell: """
        mkdir -p {out_base}/{wildcards.sample_id}/trimmed

        porechop -i {input.ONT} --format fastq.gz -t 4 -o {output.ONT}
            
        """
## The data I have seen thus far already have UMIs as part of their name:
# @A00606:487:H75CJDSX3:1:1622:10529:12493:AACCACACA 1:N:0:GATCCATG+CAACTCCA where "AACCACACA" is the UMI

rule downsample:
    input: 
        R1 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_1.fq.gz",
        R2 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_2.fq.gz"
    output: 
        R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
        R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz"
    conda: "configs/conda.yaml"
    threads: 1
    shell: """

    seqtk sample -s100 {input.R1} {sample_reads} | gzip -cvf > {output.R1}
    seqtk sample -s100 {input.R2} {sample_reads} | gzip -cvf > {output.R2}
        
    """

#Kraken2
# TODO: Should kraken use the full read set or the downsampled? ####
rule kraken2:
    input:         
        R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
        R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz",
        ONT = "{out_base}/{sample_id}/trimmed/{sample_id}_ONT_trimmed.fq.gz"
    output: 
        paired = "{out_base}/{sample_id}/{sample_id}_paired_kraken2_reads_report.txt",
        ONT = "{out_base}/{sample_id}/{sample_id}_ONT_kraken2_reads_report.txt"
    threads: 8
    conda: "configs/conda.yaml"
    shell: """

        kraken2 --db {kraken2_db} --report {output.paired} --threads 8 --paired {input.R1} {input.R2}
        kraken2 --db {kraken2_db} --report {output.ONT} --threads 8 {input.ONT}
        
    """
if config['option'] == 'UMI':

    # Assembly using unicycler
    rule assemble:
        input: 
            R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
            R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz",
            ONT = "{out_base}/{sample_id}/trimmed/{sample_id}_ONT_trimmed.fq.gz"
        output: 
            contigs = "{out_base}/{sample_id}/consensus/{sample_id}_contigs.fasta"
        conda: "configs/conda.yaml"
        threads: 8
        shell: """

            mkdir -p {out_base}/{wildcards.sample_id}/assembly
            unicycler --min_fasta_length 500 -1 {input.R1} -2 {input.R2} -l {input.ONT} -o {out_base}/{wildcards.sample_id}/assembly --threads 8

            cp {out_base}/{wildcards.sample_id}/assembly/assembly.fasta {output.contigs}      
            
            """

    # Map reads to assembly to utilise UMIs

    rule bwa_map:
        input:
            R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
            R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz",
            contigs = "{out_base}/{sample_id}/consensus/{sample_id}_contigs.fasta"
        output:
            "{out_base}/{sample_id}/mapped_reads/{sample_id}.bam"
        conda: "configs/conda.yaml"
        threads: 8
        shell: """
            bwa index {input.contigs}
            bwa mem {input.contigs} {input.R1} {input.R2} -t 8 | samtools sort > {output}
            samtools index {output}

            """

    ## The data I have seen thus far already have UMIs as part of their name:
    # @A00606:487:H75CJDSX3:1:1622:10529:12493:AACCACACA 1:N:0:GATCCATG+CAACTCCA where "AACCACACA" is the UMI        
    rule umi_tools:
        input:
            "{out_base}/{sample_id}/mapped_reads/{sample_id}.bam"
        output:
            "{out_base}/{sample_id}/mapped_reads/{sample_id}_deduplicated.bam"
        conda: "configs/conda.yaml"
        shell: """
            umi_tools dedup -I {input} --output-stats={out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}_output_stats --umi-separator=':' -S {output}
            samtools index {output}

            """

    rule consensus:
        input:
            mapping = "{out_base}/{sample_id}/mapped_reads/{sample_id}_deduplicated.bam",
            contigs = "{out_base}/{sample_id}/consensus/{sample_id}_contigs.fasta"
        output:
            "{out_base}/{sample_id}/{sample_id}_consensus.fasta"
        conda: "configs/conda.yaml"
        shell: """

            bcftools mpileup --fasta-ref {input.contigs} {input.mapping} | bcftools call -m -o {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf
            bgzip -f {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf
            tabix {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf.gz
            bcftools consensus --fasta-ref {input.contigs} {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf.gz -o {output}

            """

    rule mapping_qc:
        input:
            "{out_base}/{sample_id}/mapped_reads/{sample_id}_deduplicated.bam"
        output:
            "{out_base}/{sample_id}/qualimapReport.html"        
        conda: "configs/qc.yaml"
        threads: 4
        shell: """

            qualimap bamqc -bam {input} -nt 4 -outdir {out_base}/{wildcards.sample_id}/


            """

else:
    # Assembly using unicycler
    rule assemble:
        input: 
            R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
            R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz",
            ONT = "{out_base}/{sample_id}/trimmed/{sample_id}_ONT_trimmed.fq.gz"
        output: 
            contigs = "{out_base}/{sample_id}/{sample_id}_consensus.fasta"
        conda: "configs/conda.yaml"
        threads: 8
        shell: """

            mkdir -p {out_base}/{wildcards.sample_id}/assembly
            unicycler --min_fasta_length 500 -1 {input.R1} -2 {input.R2} -l {input.ONT} -o {out_base}/{wildcards.sample_id}/assembly --threads 8

            cp {out_base}/{wildcards.sample_id}/assembly/assembly.fasta {output.contigs}      
            
            """

rule plasmidfinder:
    input: 
        "{out_base}/{sample_id}/{sample_id}_consensus.fasta" 
    output: 
        "{out_base}/{sample_id}/plasmidfinder/{sample_id}_data.json" 
    conda: "configs/plasmidfinder.yaml"
    threads: 1
    shell: """

            mkdir -p {out_base}/{wildcards.sample_id}/plasmidfinder

            plasmidfinder.py -i {input} -o {out_base}/{wildcards.sample_id}/plasmidfinder -p {plasmidfinder_db} 
            mv {out_base}/{wildcards.sample_id}/plasmidfinder/data.json {output} 
            
            """

rule qc_assemble:
    input: 
        "{out_base}/{sample_id}/{sample_id}_consensus.fasta"
    output: 
        assembly_stats = "{out_base}/{sample_id}/assembly/{sample_id}_report.txt"
    conda: "configs/qc.yaml"
    threads: 1
    shell: """

        quast -o {out_base}/{wildcards.sample_id}/assembly/ {input}     
        cp {out_base}/{wildcards.sample_id}/assembly/report.txt {output.assembly_stats}    
        """  



rule annotate_genes:
    input:
        "{out_base}/{sample_id}/{sample_id}_consensus.fasta"
    output:
        gff = "{out_base}/{sample_id}/{sample_id}_consensus.gff",
        gbk = "{out_base}/{sample_id}/{sample_id}_consensus.gbk"
    conda: "configs/prokka.yaml"
    threads: 8
    shell: """
        mkdir -p {out_base}/{wildcards.sample_id}/prokka
        prokka --outdir {out_base}/{wildcards.sample_id}/prokka --cpu 8 --force --prefix {wildcards.sample_id} {input}
        cp {out_base}/{wildcards.sample_id}/prokka/{wildcards.sample_id}.gff {output.gff}
        cp {out_base}/{wildcards.sample_id}/prokka/{wildcards.sample_id}.gbk {output.gbk}

        """

rule multiqc:
    input:
        expand("{out_base}/{sample_id}/{sample_id}_consensus.gff", out_base = out_base, sample_id = df["sample_id"])
    output:
        "{out_base}/multiqc_report.html"
    conda: "configs/qc.yaml"
    threads: 1
    shell: """
        multiqc -d {out_base} -o {out_base}

        """

