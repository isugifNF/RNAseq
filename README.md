# RNAseq

RNAseq from fastq reads to gene counts/quantification

```
nextflow run main.nf --help
```

```
N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [curious_curran] - revision: 7d9f51e73b
Usage:
   The typical command for running the RNAseq pipeline is as follows:
   nextflow run main.nf --reads "*_{R1,R2}.fastq.gz" --genome GENOME.fasta --genome_gff GENOME.gff --genome_cdna CDNA.fasta  -profile singularity

   Mandatory arguments:
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
    --genome                Reference genome fasta file
    --genome_gff            Reference gff file
    --genome_cdna           Reference cdna file

   Optional configuration arguments:
    --methods               Select one or more RNAseq counting methods [default: 'gsnap,hisat2,kallisto,salmon']
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
mkdir 00_Raw-Data; cd 00_Raw-Data;
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
cd 00_Raw-Data

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

<details><summary>See example run on <b>Ceres HPC</b> - last update 24 June 2021 </summary>

Ran in 1 hour and 48 minutes. The below commands were wrapped in a slurm script.

```
# Ceres Modules
module load nextflow
module load parallel 
module load fastqc
module load python_3
module load salmon/0.10.1
module load kallisto/0.42.4
module load hisat2/2.2.0
module load gmap_gsnap/2020-04-08
module load samtools
module load subread/2.0.2

# Main command
nextflow run main.nf \
  --reads "00_Raw-Data/*{1,2}.fastq.gz" \
  --genome "00_Raw-Data/*_genomic.fna.gz" \
  --genome_gff "00_Raw-Data/*.gff.gz" \
  --genome_cdna "00_Raw-Data/*_rna.fna.gz" \
  --queueSize 25 \
  -profile slurm \
  -resume
```

Progress messages

```
N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [elated_kilby] - revision: 7d9f51e73b
executor >  slurm (117)
[5e/2c59b1] process > fastqc (batched)               [100%] 5 of 5 ✔
[90/ed372f] process > multiqc                        [100%] 1 of 1 ✔
[b7/06aa08] process > kallisto_index (GCF_9021671... [100%] 1 of 1 ✔
[62/2027e7] process > kallisto_quant (SRR1573520)    [100%] 18 of 18 ✔
[3e/ebd42f] process > salmon_index (GCF_902167145... [100%] 1 of 1 ✔
[d9/03b7df] process > salmon_quant (SRR1573520)      [100%] 18 of 18 ✔
[4c/ec81ca] process > gsnap_index (GCF_902167145)    [100%] 1 of 1 ✔
[36/237247] process > gsnap_align (SRR1573520)       [100%] 18 of 18 ✔
[5f/54c2e9] process > featureCounts_gene (SRR1573... [100%] 18 of 18 ✔
[c7/f2c642] process > featureCounts_mRNA (SRR1573... [100%] 18 of 18 ✔
[db/3e6c52] process > featureCounts_geneMult (SRR... [100%] 18 of 18 ✔
Completed at: 24-Jun-2021 13:07:28
Duration    : 1h 48m 38s
CPU hours   : 11.5
Succeeded   : 117
```

</details>