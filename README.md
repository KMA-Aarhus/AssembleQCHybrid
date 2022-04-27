# AssembleQCHybrid
This workflow uses output from both paired-end and nanopore sequencing to produce a high quality assembly.
Input: 
- Raw nanopore fastq files as output from MinKnow. Raw nanopore fastq from "fastq_pass" must be joined. Example: 
```
find /fastq_pass/barcode01 -name '\*fastq.gz' -exec cat {} + > barcode01.fastq.gz \
```
Output:  
- QC report. 
- Taxonomic analysis. 
- Genome Assembly. 
- Genome annotation. 
## Installation. 
Requires conda, openpyxl and snakemake to run.  
https://docs.conda.io/en/latest/miniconda.html. 
https://openpyxl.readthedocs.io/en/stable/.  
https://snakemake.readthedocs.io/en/stable/. 

### Install with git
```
git clone https://github.com/KMA-Aarhus/AssembleQCHybrid.git
```
## DAG of workflow
![assembleqchybrid](https://user-images.githubusercontent.com/90172976/165513150-c4cd306d-ab01-4964-b82c-1b0dfd637550.png)

## How to run
Set up an alias for snakemake to run on slurm:
```
alias snakeslurm='mkdir -p logs/old; mv logs/*.{err,out} logs/old 2> /dev/null; snakemake --profile configs/slurm --use-conda --conda-frontend mamba'
```
Navigate to AssemblyQCHybrid directory where the snakefile is.  


### To run 
```
snakeslurm --config R1="path_to_/R1.fastq" R2="path_to_/R2.fastq" ONT="path_nanopore.fastq"
Optionally, an out_base can be provided. as out_base="output". The workflow will then save output files to this location. Default="output".
If your paired-end reads contain UMIs, specify this using option="UMI"

Additional optional options:
kraken2_db # Location of the kraken2 database. Default="/project/ClinicalMicrobio/faststorage/database/kraken2_20210423"
```
### Output
Output can be found in the specified output directory. Output contains:
* A multi-qc report for all samples in the run. The report contains read statistics, assembly statistics, taxonomic analysis, genome coverage and more.
* One folder per sample.
  * "sample_name"_consensus.fasta: assembled and polished genome
  * "sample_name"_consensus.gff: genome annotations
  * "sample_name"_consensus.gbk: annotated genome in genbank format
  * "sample_name"_kraken2_reads_report.txt: the output from kraken2. The same information can be found in the multi-qc report.
  * Addtional output for debugging can be found in:
    * qc
    * plasmidfinder
    * prokka
    * assembly
    * trimmed
    * sampled

