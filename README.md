# RNAseq

RNAseq from fastq reads to gene counts/quantification

```
nextflow run main.nf --help
```

```
N E X T F L O W  ~  version 21.04.0
Launching `main.nf` [berserk_celsius] - revision: 8a7622a6be

   Usage:
   The typical command for running the RNAseq pipeline is as follows:
   nextflow run main.nf --reads "*_{R1,R2}.fastq.gz" --genome GENOME.fasta --genome_gff GENOME.gff --genome_cdna CDNA.fasta  -profile singularity

   Mandatory arguments:
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
    --genome                Reference genome fasta file
    --genome_gff            Reference gff file
    --genome_cdna           Reference cdna file

   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, slurm [default:local]
    --fastqc_app            Link to fastqc executable [default: 'fastqc']
    --multiqc_app           Link to multiqc executable [default: 'multiqc']
    --parallel_app          Link to parallel executable [default: 'parallel']
    --kallisto_app          Link to kallisto executable [default: 'kallisto']
    --salmon_app            Link to salmon executable [default: 'salmon']
    --gsnap_app             Link to gsnap executable [default: 'gsnap']
    --featureCounts_app     Link to featureCounts executable [default: 'featureCounts']

   Optional other arguments:
    --threads               Threads per process [default:4 for local, 16 for slurm]
    --queueSize             Maximum jobs to submit to slurm [default:20]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help
```
## Example Dataset

<details><summary>Fetch Maize Reads</summary>

```
mkdir 00_INPUT; cd 00_INPUT;
# ==== Fetch Reads
# === Blade Tissue
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573504/SRR1573504_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573505/SRR1573505_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573506/SRR1573506_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573506/SRR1573506_2.fastq.gz

# === Blade mutant L1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573519/SRR1573519_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573519/SRR1573519_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573520/SRR1573520_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573520/SRR1573520_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573521/SRR1573521_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573521/SRR1573521_2.fastq.gz

# === Ligule Tissue
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573507/SRR1573507_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573507/SRR1573507_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573508/SRR1573508_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/008/SRR1573508/SRR1573508_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573509/SRR1573509_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/009/SRR1573509/SRR1573509_2.fastq.gz

# === Ligule mutant L1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573522/SRR1573522_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573522/SRR1573522_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573523/SRR1573523_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/003/SRR1573523/SRR1573523_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573524/SRR1573524_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/004/SRR1573524/SRR1573524_2.fastq.gz

# === Sheath Tissue
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573510/SRR1573510_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/000/SRR1573510/SRR1573510_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573511/SRR1573511_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/001/SRR1573511/SRR1573511_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573512/SRR1573512_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/002/SRR1573512/SRR1573512_2.fastq.gz

# === Sheath mutant L1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573525/SRR1573525_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR1573525/SRR1573525_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573526/SRR1573526_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/006/SRR1573526/SRR1573526_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573527/SRR1573527_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/007/SRR1573527/SRR1573527_2.fastq.gz
```

</details>

Reference files

Fetch `Zea mays` reference files from NCBI

* https://www.ncbi.nlm.nih.gov/assembly/GCF_902167145.1/

```
mkdir 02_Genome; cd 02_Genome
 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_rna.fna.gz
```

## Example Run

```
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [naughty_keller] - revision: 6370a27351
executor >  slurm (154)
[65/9f5b8d] process > fastqc (batched)               [100%] 6 of 6 ✔
[a6/f76504] process > multiqc                        [100%] 1 of 1 ✔
[3f/cb5455] process > kallisto_index (GCF_9021671... [100%] 1 of 1 ✔
[96/940514] process > kallisto_quant (SRR1573504)    [100%] 24 of 24 ✔
[2d/ab4e35] process > salmon_index (GCF_902167145... [100%] 1 of 1 ✔
[b5/260d0a] process > salmon_quant (SRR1573509)      [100%] 24 of 24 ✔
[42/d3600c] process > gsnap_index (GCF_902167145)    [100%] 1 of 1 ✔
[17/b68c59] process > gsnap_quant (SRR1573511)       [100%] 24 of 24 ✔
[51/de2e96] process > featureCounts_gene (SRR1573... [100%] 24 of 24 ✔
[5b/3c04f4] process > featureCounts_mRNA (SRR1573... [100%] 24 of 24 ✔
[7c/7ccc62] process > featureCounts_geneMult (SRR... [100%] 24 of 24 ✔
Completed at: 22-Mar-2021 04:51:51
Duration    : 4h 20m 7s
CPU hours   : 71.3
Succeeded   : 154
```

[timeline.html](https://isugifnf.github.io/RNAseq/timeline.html) | [report.html](https://isugifnf.github.io/RNAseq/report.html)

```
RNASeq_Results/
  |_ 01_Quality-Control/     # fastqc and multiqc report
  |_ 03_GSNAP                # gene counts from GSNAP
  |_ 03_Kallisto             # transcript counts from kallisto (ready for sleuth analysis)
  |_ 03_Salmon               # transcript counts from salmon (ready for DESeq2 analysis)
  |_ report.html             # memory usage report
  |_ timeline.html           # runtime report
```
