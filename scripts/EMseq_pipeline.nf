#!/usr/bin/env nextflow

// Define parameters
params.highreads = "../data/high_yield/*_{R1,R2}*.fastq.gz"
params.lowreads = "../data/low_yield/*_{R1,R2}*.fastq.gz"
params.hightrimmed = "../output/hightrimmedfastq"
params.hightrimmedreads = "../output/hightrimmedfastq/*_{R1,R2}*.fq.gz"
params.lowtrimmed = "../output/lowtrimmedfastq"
params.genomeDir = "../data/genome_test"
params.genomeHighBAM = "../output/high_bam"
params.genomeLowBAM = "../output/low_bam"

// Run trim galore process
process hightrimGalore {
	// Copy the output to the hightrimmedfastq directory from work directory
	publishDir("${params.hightrimmed}", mode: "copy")
	// Define the number of CPUs to use
	cpus 4
	// Define the memory to use
	memory 6.GB
	// Define the input
	input:
	tuple val(pair_id), file(fastq)

	output:
	path "*"

	script:
	"""
	trim_galore --paired --clip_R1 2 --clip_R2 2 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --fastqc --quality 20 --polyA --illumina --length 100 $fastq
	"""
}

process lowtrimGalore {
	// Copy the output to the lowtrimmedfastq directory from work directory
	publishDir("${params.lowtrimmed}", mode: "copy")
	// Define the number of CPUs to use
	cpus 4
	// Define the memory to use
	memory 6.GB
	// Define the input
	input:
	tuple val(pair_id), file(fastq)

	output:
	path "*"

	script:
	"""
	trim_galore --paired --clip_R1 2 --clip_R2 2 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --fastqc --quality 20 --polyA --illumina --length 100 $fastq
	"""
}

// Run Bismark genome preparation process
process genomePreparation {
	// Copy the output to the genome directory from work directory
	publishDir("${params.genomeDir}", mode: 'copy')
	// Define the number of CPUs to use
	cpus 4
	// Define the memory to use
	memory 6.GB
	// Define the input
	input:
	path genome

	script:
	"""
	bismark_genome_preparation $genome
	"""
}

// Run Bismark high reads alignment process
process highRead_alignment {
	// Copy the output to the highReads_bam directory from work directory
	publishDir("${params.genomeHighBAM}", mode: 'copy')
	// Define the number of CPUs to use
	cpus 4
	// Define the memory to use
	memory 6.GB
	// Define the input
	input:
	tuple val(read_id), path(fastq1), path(fastq2)
	// genome directory
	path genomeDir

	output:
	path "*_bismark_bt2_pe.bam", emit: aligned_bam_files
	path "*_bismark_bt2_PE_report.txt", emit: bismark_reports

	script:
	"""
	bismark --un --ambiguous --genome $genomeDir -1 ${fastq1} -2 ${fastq2}
	"""
}

workflow {
	// Define channels
	high_trim = Channel.fromFilePairs(params.highreads)
	low_trim = Channel.fromFilePairs(params.lowreads)
	genome_prep = Channel.fromPath(params.genomeDir)
	high_trimmed = Channel.fromFilePairs(params.hightrimmedreads, flat: true)

	// run trim galore processes
	hightrimGalore(high_trim)
	lowtrimGalore(low_trim)

	// run genome preparation process
	genomePreparation(genome_prep)

	// run read alignment process for high reads
	highRead_alignment(high_trimmed , genome_prep)
}