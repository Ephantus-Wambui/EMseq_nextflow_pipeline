# EMseq nextflow pipeline

This pipeline is designed to process raw data from EM-seq experiments. It is based on the [nf-core](https://nf-co.re/) framework.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

Ensure also on your system you have installed fastqc, multiqc, trim galore and bismark.

2. Install [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
3. Install [`multiqc`](https://multiqc.info/)
4. Install [`trim galore`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
5. Install [`bismark`](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)

## Clone repo

```bash
git clone [`EMseq pipeline`](https://github.com/Ephantus-Wambui/EMseq_nextflow_pipeline.git)
```

## Run pipeline

Before running the pipeline ensure you have the following files in the `data` directory:

1. genome_test directory which contains TMEB117_chr16.fasta reference genome

2. high yield and low yield fastq files

After ensuring everything is in place, activate the conda environment which contains fastqc, multiqc, trim galore and bismark dependencies.

```bash
conda activate EMseq # activate conda environment
```

Before running the EMseq_pipeline.nf script, first run the QC pipeline to check the quality of the fastq file.

** Note: cd into scripts folder **

```bash
nextflow run EMseq_fastqc.nf
```

After running the QC pipeline, run the EMseq pipeline to align the reads to the reference genome and generate methylation calls.

```bash
nextflow run EMseq_pipeline.nf
```
** Note: Adjust trim galore parameters according to the fastqc results and then run the pipeline. **

## Output

The pipeline will generate the following files directories in the output directory:

1. Both high yield and low yield directories which will contain individual fastqc reports of the fastq files.

2. Both high yield and low yield directories which will contain multiqc reports of the fastqc reports.

3. Both high yield and low yield directories which will contain trimmed fastq files.

4. Both high yield and low yield directories which will contain bismark alignment reports.

5. Both high yield and low yield directories which will contain bismark methylation calls.

## Pipeline summary

1. QC pipeline: This pipeline will generate fastqc reports of the fastq files and a multiqc report of the fastqc reports.

2. EMseq pipeline: This pipeline will align the reads to the reference genome and generate methylation calls.

## Contributors

1. [Ephantus Wambui](https://github.com/Ephantus-Wambui/)
