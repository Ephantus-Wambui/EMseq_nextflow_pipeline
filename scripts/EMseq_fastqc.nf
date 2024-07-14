#!/usr/bin/env nextflow

params.highreads = "../data/high_yield/*.fastq.gz"

params.lowreads = "../data/low_yield/*.fastq.gz"

params.fastqc_highreports = "../output/fastqc_highresults"

params.fastqc_lowreports = "../output/fastqc_lowresults"

params.multiqc_highreports = "../output/multiqc_highresults"

params.multiqc_lowreports = "../output/multiqc_lowresults"

// Run FastQC process
process highfastqc {
    publishDir("${params.fastqc_highreports}", mode: "copy")
    input:
    path fastq

    output:
    path "*"

    script:
    """
    fastqc $fastq
    """
}

process lowfastqc {
    publishDir("${params.fastqc_lowreports}", mode: "copy")
    input:
    path fastq
    output:
    path "*"
    script:
    """
    fastqc $fastq
    """
}

// Run MultiQC process
process highmultiqc {
    publishDir("${params.multiqc_highreports}", mode: "copy")
    input:
    path fastqc_highresults

    output:
    path "*"

    script:
    """
    multiqc .
    """
}

process lowmultiqc {
    publishDir("${params.multiqc_lowreports}", mode: "copy")
    input:
    path fastqc_lowreports
    output:
    path  "*"
    script:
    """
    multiqc .
    """
}

workflow {
	high_fastq = Channel.fromPath(params.highreads)
	low_fastq = Channel.fromPath(params.lowreads)

	highfastqc(high_fastq) | collect | highmultiqc
	lowfastqc(low_fastq) | collect | lowmultiqc
}
