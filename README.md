# RSVartic_wrapper pipeline

A Snakemake-based pipeline for automating the detection of RSV subtypes, genome assembly using reference sequences, and mutation analysis from fragmented rapid barcoding nanopore sequencing data.


## Overview
This workflow integrates several powerful tools:

BLASTn for identifying RSV subtypes.
ARTIC v1.5 from the Field Bioinformatics toolkit for genome assembly and polishing.
Custom in-house methods for detecting severe mutations, based on literatures.

## Dependencies

Ensure the following tools and libraries are installed:
- seqkit
- blast=2.16
- artic1.5
- python
- biopython=1.84
- pysam

## Installation and Setup

Follow these steps to clone and install the repository (Linux-based systems):

	```
	git clone ssh https://git@github.com:aradahir/RSVartic_wrapper.git
	cd .\RSVartic_wrapper\
	conda env create -f environment.yml

	```

Activate the snakemake environment for running the pipeline.
	
	```
	conda activate artic1.5
	
	```

## Usage

 - Place raw nanopore sequencing fastq folder in ./data/raw/ directory.
 	 
 - Open the environment
 
 		```
		conda activate snakemake
		```
 - Run the pipeline
 		-
 		thread can be adjusted by changing from 1 into the specific number of thread

 		```
 		snakemake -j 1
		```

# result interpretation

there are 4 output in data folder from this pipeline
1. concatenated_files
	- the concatenate of fastq file from raw data named as the given name from samplesheet.csv
2. fasta_blastn
	- the fasta file from each contig
	- the .txt file result of each contig searching agiast blastn
3. subtype folder
	- the RSVA or RSVB folder that contain the concatenated fastq files.

there are 7 output in results folder from this pipeline 
1. data: 
	- subtypes folder including RSVA, RSVB, and unknown which insides contain the concatenate fastq file
2. artic:
	- all results from artic workflow.
3. depth:
	- depth coverage plot before and after using artic pipeline.
4. F_protiens:
	- a folder containing only amino acid sequences of F protiens from each sequence.
	- a folder of mutations derived from comparing the sample against literatures.
	

[!WARNING]
This pipeline developed for nanopore sequencing prepared by fragmentation of rapid barcoding library preparation protocol. 
