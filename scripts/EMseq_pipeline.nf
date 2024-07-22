#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.highreads = "../data/high_yield/*_{R1,R2}*.fastq.gz"
params.lowreads = "../data/low_yield/*_{R1,R2}*.fastq.gz"
params.hightrimmed = "../output/hightrimmedfastq"
params.lowtrimmed = "../output/lowtrimmedfastq"
params.genomeDir = "../data/genome_test"
params.genomeHighBAM = "../output/high_bam"
params.genomeLowBAM = "../output/low_bam"
params.highDedups = "../output/highDedups"
params.lowDedups = "../output/lowDedups"
params.highMeth = "../output/highMeth"
params.lowMeth = "../output/lowMeth"

// Define processes
process hightrimGalore {
    publishDir(params.hightrimmed, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    tuple val(pair_id), file(fastq)

    output:
    tuple val(pair_id), path("*.fq.gz"), emit:high_trim_galore

    script:
    """
    trim_galore --paired --clip_R1 2 --clip_R2 2 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --fastqc --quality 20 --polyA --illumina --length 100 $fastq
    """
}

process lowtrimGalore {
    publishDir(params.lowtrimmed, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    tuple val(pair_id), file(fastq)

    output:
    tuple val(pair_id), path("*.fq.gz"), emit:low_trim_galore

    script:
    """
    trim_galore --paired --clip_R1 2 --clip_R2 2 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --fastqc --quality 20 --polyA --illumina --length 100 $fastq
    """
}

process genomePreparation {
    publishDir(params.genomeDir, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    path genome

    script:
    """
    bismark_genome_preparation $genome
    """
}

process highRead_alignment {
    publishDir(params.genomeHighBAM, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    tuple val(pair_id), path(fastq)
    path genome

    output:
    path "*_bismark_bt2_pe.bam", emit:high_bam
    path "*_bismark_bt2_PE_report.txt", emit: bismark_reports

    script:
    """
    bismark --un --ambiguous --genome $genome -1 ${fastq[0]} -2 ${fastq[1]}
    """
}

process lowRead_alignment {
    publishDir(params.genomeLowBAM, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    tuple val(pair_id), path(fastq)
    path genome

    output:
    path "*_bismark_bt2_pe.bam", emit:low_bam
    path "*_bismark_bt2_PE_report.txt", emit: bismark_reports

    script:
    """
    bismark --un --ambiguous --genome $genome -1 ${fastq[0]} -2 ${fastq[1]}
    """
}

process highDedups {
    publishDir(params.highDedups, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    path bam

    output:
    path "*_bismark_bt2*.deduplicated.bam", emit:high_dedups
    path "*deduplication_report.txt", emit: deduplication_reports

    script:
    """
    deduplicate_bismark $bam
    """
}

process lowDedups {
    publishDir(params.lowDedups, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    path bam

    output:
    path "*_bismark_bt2*.deduplicated.bam", emit:low_dedups
    path "*deduplication_report.txt", emit: deduplication_reports

    script:
    """
    deduplicate_bismark $bam
    """
}

process highMethExtraction {
    publishDir(params.highMeth, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    path bam

    output:
    path("*")

    script:
    """
    bismark_methylation_extractor $bam
    """
}

process lowMethExtraction {
    publishDir(params.lowMeth, mode: 'copy')
    cpus 4
    memory 6.GB

    input:
    path bam

    output:
    path("*")

    script:
    """
    bismark_methylation_extractor $bam
    """
}

workflow {
    // Define channels
    high_trim = Channel.fromFilePairs(params.highreads, size: 2)
    low_trim = Channel.fromFilePairs(params.lowreads, size: 2)
    genome_prep = Channel.fromPath(params.genomeDir)

    // Run genome preparation process once
    genomePreparation(genome_prep)

    // Run trim galore processes
    hightrimGalore(high_trim)
    lowtrimGalore(low_trim)

    // Run read alignment process for high and low reads
    highRead_alignment(hightrimGalore.out.high_trim_galore, genome_prep)
    lowRead_alignment(lowtrimGalore.out.low_trim_galore, genome_prep)

    // Run deduplication process for high and low reads
    highDedups(highRead_alignment.out.high_bam)
    lowDedups(lowRead_alignment.out.low_bam)

    // Run methylation extraction process for high and low reads
    highMethExtraction(highDedups.out.high_dedups)
    lowMethExtraction(lowDedups.out.low_dedups)
}
