#!/usr/bin/env nextflow

// Define parameters
params.highreads = "../data/high_yield/*_{R1,R2}*.fastq.gz"
params.lowreads = "../data/low_yield/*_{R1,R2}*.fastq.gz"
params.hightrimmed = "../output/hightrimmedfastq"
params.hightrimmedreads = "../output/hightrimmedfastq/*_{R1,R2}*.fq.gz"
params.lowtrimmedreads = "../output/lowtrimmedfastq/*_{R1,R2}*.fq.gz"
params.genomeDir = "../data/genome_test"
params.genomeHighBAM = "../output/high_bam"
params.genomeLowBAM = "../output/low_bam"
params.highBamFile = "../output/high_bam/*.bam"
params.lowBamFile = "../output/low_bam/*.bam"
params.highDedups = "../output/highDedups"
params.lowDedups = "../output/lowDedups"
params.highDedupsFile = "../output/highDedups/*.bam"
params.lowDedupsFile = "../output/lowDedups/*.bam"
params.highMeth = "../output/highMeth"
params.lowMeth = "../output/lowMeth"

// Create trim galore process
process hightrimGalore {
	// Copy the output to the hightrimmedfastq directory from work directory
	publishDir("${params.hightrimmed}", mode: "copy")
	// Define the number of CPUs to use
	cpus 2
	// Define the memory to use
	memory 3.GB
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
	cpus 2
	// Define the memory to use
	memory 3.GB
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

// Create Bismark genome preparation process
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

// Create Bismark high reads alignment process
process highRead_alignment {
	// Copy the output to the highReads_bam directory from work directory
	publishDir("${params.genomeHighBAM}", mode: 'copy')
	// Define the number of CPUs to use
	cpus 2
	// Define the memory to use
	memory 3.GB
	// Define the input
	input:
	tuple val(pair_id), file(fastq)
	// Genome directory
	path genomeDir

	output:
	path "*"

	script:
	"""
	bismark --un --ambiguous --genome $genomeDir -1 ${fastq[0]} -2 ${fastq[1]}
	"""
}

// Create Bismark low reads alignment process
process lowRead_alignment {
    // Copy the output to the highReads_bam directory from work directory
    publishDir("${params.genomeLowBAM}", mode: 'copy')
    // Define the number of CPUs to use
    cpus 2
    // Define the memory to use
    memory 3.GB
    // Define the input
    input:
    tuple val(pair_id), file(fastq)
	// Genome directory
	path genomeDir

    output:
    path "*"

    script:
    """
    bismark --un --ambiguous --genome $genomeDir -1 ${fastq[0]} -2 ${fastq[1]}
    """
}

// Create deduplication process
process highDedups {
	// Copy output to the highDedups directory from work directory
	publishDir("${params.highDedups}", mode: 'copy')
	// Define the number of CPUs to use
	cpus 2
	// Define the memory to use
	memory 3.GB
	// Define the input
	input:
	path bam

	output:
	path "*"

	script:
	"""
	deduplicate_bismark $bam
	"""
}

// Create low deduplication process
process lowDedups {
	// Copy output to the lowDedups directory from work directory
	publishDir("${params.lowDedups}", mode: 'copy')
	// Define the number of CPUs to use
	cpus 2
	// Define the memory to use
	memory 3.GB
	// Define the input
	input:
	path bam

	output:
	path "*"

	script:
	"""
	deduplicate_bismark $bam
	"""
}

// Create high methylation extraction process
process highMethExtraction {
	// Copy output to the highMeth directory from work directory
	publishDir("${params.highMeth}", mode: 'copy')
	// Define the number of CPUs to use
	cpus 2
	// Define the memory to use
	memory 3.GB
	// Define the input
	input:
	path bam

	output:
	path "*"

	script:
	"""
	bismark_methylation_extractor $bam
	"""
}

// Create low methylation extraction process
process lowMethExtraction {
	// Copy output to the lowMeth directory from work directory
	publishDir("${params.lowMeth}", mode: 'copy')
	// Define the number of CPUs to use
	cpus 2
	// Define the memory to use
	memory 3.GB
	// Define the input
	input:
	path bam

	output:
	path "*"

	script:
	"""
	bismark_methylation_extractor $bam
	"""
}

workflow {
	// Define channels
	high_trim = Channel.fromFilePairs(params.highreads)
	low_trim = Channel.fromFilePairs(params.lowreads)
	genome_prep = Channel.fromPath(params.genomeDir)
	high_trimmed = Channel.fromFilePairs(params.hightrimmedreads)
	low_trimmed = Channel.fromFilePairs(params.lowtrimmedreads)
	high_bam = Channel.fromPath(params.highBamFile)
	low_bam = Channel.fromPath(params.lowBamFile)
	high_dedups = Channel.fromPath(params.highDedupsFile)
	low_dedups = Channel.fromPath(params.lowDedupsFile)

	// Run trim galore processes
	// hightrimGalore(high_trim)
	// lowtrimGalore(low_trim)

	// Run genome preparation process
	// genomePreparation(genome_prep)

	// Run read alignment process for high reads
	// highRead_alignment(high_trimmed, genome_prep)
	// lowRead_alignment(low_trimmed, genome_prep)

	// Run deduplication process
	// highDedups(high_bam)
	// lowDedups(low_bam)

	// Run methylation extraction process
	highMethExtraction(high_dedups)
	lowMethExtraction(low_dedups)
}
