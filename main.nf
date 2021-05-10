#! /usr/bin/env nextflow
/* Name: RNASeq pipeline
 * Auth: Jennifer Chang
 * Date: 2021/03/14
 * Desc: Given a reference genome and paired reads, create RNAseq counts files
 */

nextflow.enable.dsl=2

// === Define Processes

process fastqc {
    tag "batched"
    label 'fastqc'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'
//    module 'fastqc:parallel'

    publishDir "${params.outdir}/01_Quality-Control", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

    input:
    path(read)

    output:
    path("*_fastqc.{zip,html}")

    script:
    """
    #! /usr/bin/env bash
    $parallel_app -j8 "$fastqc_app {1}" ::: $read
    """
}

process multiqc {
    label 'multiqc'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/01_Quality-Control", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

    input:
    path fastqc_htmls

    output:
    path "multiqc_report.html"

    script:
    """
    #! /usr/bin/env bash
    $multiqc_app .
    """
}

process kallisto_index {
    tag "$genome_cdna"
    label 'kallisto'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_Kallisto", mode: 'copy'

    input:
    path(genome_cdna)

    output:
    path("${genome_cdna.simpleName}.idx")

    script:
    """
    #! /usr/bin/env bash
    $kallisto_app index \
      -i ${genome_cdna.simpleName}.idx ${genome_cdna}
    """
}

process kallisto_quant {
    tag "${readname}"
    label 'kallisto'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_Kallisto", mode: 'copy'

    input:
    tuple path(genome_index), val(readname), path(read_pairs)

    output:
    path("*")

    script:
    """
    #! /usr/bin/env bash
    $kallisto_app quant \
     -i ${genome_index} \
     -o ${readname}_quant \
     -b 20 -t 16 \
     ${read_pairs}
    """
}

process salmon_index {
    tag "$genome_cdna"
    label 'salmon'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_Salmon", mode: 'copy'

    input:
    path(genome_cdna)

    output:
    tuple val("${genome_cdna.simpleName}"), path("*")

    script:
    """
    #! /usr/bin/env bash
    $salmon_app index \
      -i ${genome_cdna.simpleName} \
      -t ${genome_cdna}
    """
}

process salmon_quant {
    tag "${readname}"
    label 'salmon'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 04:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_Salmon", mode: 'copy'

    input:
    tuple val(genome_index), path(genome_index_files), val(readname), path(read_pairs)

    output:
    path("*")

    script:
    """
    #! /usr/bin/env bash
    $salmon_app quant \
     -l A -p 16 \
     --validateMappings \
     -i ${genome_index} \
     -1 ${read_pairs.get(0)} \
     -2 ${read_pairs.get(1)} \
     -o ${readname}_quant
    """
}

process gsnap_index {
    tag "${genome_gz.simpleName}"
    label 'gsnap'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_GSNAP", mode: 'copy'

    input:
    path(genome_gz)

    output:
    tuple val("${genome_gz.simpleName}"), path("*")

    script:
    """
    #! /usr/bin/env bash
    gmap_build \
     --gunzip \
     -d ${genome_gz.simpleName} \
     -D gmapdb \
     ${genome_gz}
    """
}

process gsnap_align {
    tag "${readname}"
    label 'gsnap'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 04:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_GSNAP", mode: 'copy'

    input:
    tuple val(genome_name), path(gmap_dir),val(readname), path(read_pairs)

    output:
    path("*")

    script:
    """
    #! /usr/bin/env bash
    $gsnap_app \
     --gunzip \
     -d ${genome_name} \
     -D gmapdb/ \
     -N 1 -t 16 -B 4 -m 5 \
     --input-buffer-size=1000000 \
     --output-buffer-size=1000000 \
     -A sam \
     ${read_pairs} |
     samtools view --threads 16 -bS - > ${readname.simpleName}.bam
    """
}

process featureCounts_gene {
    tag "${read_bam.simpleName}"
    label 'gsnap'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 04:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_GSNAP", mode: 'copy'

    input:
    tuple path(genome_gff), path(read_bam)

    output:
    path("*")

    script:
    """
    #! /usr/bin/env bash
    $featureCounts_app \
      -T 16 \
      -t gene \
      -g ID \
      -a ${genome_gff} \
      -o ${read_bam.simpleName}_genecounts.txt \
      ${read_bam}
    """
}

process featureCounts_mRNA {
    tag "${read_bam.simpleName}"
    label 'gsnap'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 04:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_GSNAP", mode: 'copy'

    input:
    tuple path(genome_gff), path(read_bam)

    output:
    path("*")

    script:
    """
    #! /usr/bin/env bash
    $featureCounts_app \
      -T 16 \
      -t mRNA \
      -g ID \
      -a ${genome_gff} \
      -o ${read_bam.simpleName}_mRNAcounts.txt \
      ${read_bam}
    """
}

process featureCounts_geneMult {
    tag "${read_bam.simpleName}"
    label 'gsnap'
//    executor 'slurm'
//    clusterOptions '-N 1 -n 16 -t 04:00:00  --account=isu_gif_vrsc'

    publishDir "${params.outdir}/03_GSNAP", mode: 'copy'

    input:
    tuple path(genome_gff), path(read_bam)

    output:
    path("*")

    script:
    """
    #! /usr/bin/env bash
    $featureCounts_app \
      -T 16 -M \
      -t gene \
      -g ID \
      -a ${genome_gff} \
      -o ${read_bam.simpleName}_geneMultcounts.txt \
      ${read_bam}
    """
}

// === Main Workflow
workflow {

  /* Read, reference, and gff channels */
  readsflat_ch = channel.fromPath(params.reads, checkIfExists:true)
  reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
  cdna_ch = channel.fromPath(params.genome_cdna, checkIfExists:true)
  genome_ch = channel.fromPath(params.genome, checkIfExists:true)
  gff_ch = channel.fromPath(params.genome_gff, checkIfExists:true)

  /* 01_Quality-Control */
  readsflat_ch.collate(8) | fastqc | collect | multiqc

  /* 03_Kallisto */
  cdna_ch | kallisto_index | combine(reads_ch) | kallisto_quant

  /* 03_Salmon */
  cdna_ch | salmon_index | combine(reads_ch) | salmon_quant

  /* 03_GSNAP */
  genome_ch | gsnap_index | combine(reads_ch) | gsnap_align
  gff_ch.combine(gsnap_align.out) | ( featureCounts_gene & featureCounts_mRNA & featureCounts_geneMult)
}
