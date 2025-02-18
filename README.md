# DRAGEN-GATK Workflow

## Introduction
This workflow integrates **DRAGEN** and **GATK** to process genomic sequencing data, including read alignment, base quality recalibration, and variant calling. It is optimized for seamless execution through a user interface (UI) rather than the command line.

## Features
- **UI-Driven Execution**: All inputs and parameters are configured through a graphical interface.
- **Automated Processing**: Supports alignment, recalibration, and variant calling without manual intervention.
- **Reference Genome Support**: Accepts FASTA files, including `.fa`, `.fasta`, `.fa.gz`, and `.fasta.gz` formats.
- **Sample Sheet Input**: Uses structured CSV input for managing sample metadata.
- **Variant Recalibration**: Optional processing using known variant sites in VCF format.

## Input Requirements
- **Sample Sheet**: A CSV file defining sample IDs and file paths.
- **Reference Genome**: A FASTA file, optionally compressed.
- **Known Sites VCF** (optional): Used for base quality recalibration.

## Output Overview
- **Aligned Reads**: BAM files of mapped reads.
- **Recalibrated Data**: BAM files post-recalibration (if applicable).
- **Variant Calls**: VCF files containing detected variants.

For detailed usage instructions, see [usage.md](usage.md).
