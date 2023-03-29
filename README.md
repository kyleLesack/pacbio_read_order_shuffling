# Overview

These scripts were used to evaluate the impact of FASTQ file read order on structural variation (SV) calling from long-read DNA sequencing data.

# Requirements

Note: The ./yaml directory include conda yaml environment files. I recommend using conda to install the required software using these files. The required tools and the versions we tested are listed below:

* Snakemake (v.)
* [Picard Tools](https://broadinstitute.github.io/picard/) (v.)
* SAMtools (v1.9)
	* Note: SAMTools >= v1.10 didn't work with the NGMLR alignments
* [seq-shuff](https://github.com/thackl/seq-scripts/blob/master/bin/seq-shuf)
* Aligners
	* [Minimap2](https://github.com/lh3/minimap2) (v.2.17)
	* [NGMLR](https://github.com/philres/ngmlr) (v.0.2.7)
	* [pbmm2](https://github.com/PacificBiosciences/pbmm2) (v.1.7.0)
* SV Callers
	* [pbsv](https://github.com/PacificBiosciences/pbsv) (v2.8.0)
	* [Sniffles](https://github.com/fritzsedlazeck/Sniffles) (v2.0.6)
	* [SVIM](https://github.com/eldariont/svim) (v2.0.0)
* *Caenorhabditis elegans* reference genome and associated files
 	* FASTQ index
	* sawriter suffix array files and NGMLR indexes are required for pbmm2/pbsv

# Snakemake pipeline

The ./Snakefile contains the code to run the pipeline. The pipeline includes several computational steps that were performed on a high performance computing cluster. The default resource requirements that we provided to our SLURM scheduler are specified in the ./slurm/config.yaml file. The default values were overwritten for certain steps, which are specified in certain Snakefile rules under threads and resources.

# Data Not Included Here

* The sawriter files were not included due to size
	* 0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa-enc.2.ngm
	* 0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa-ht-13-2.2.ngm
	* 0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa.sa


# Test Data

We used PacBio sequencing data from the [Caenorhabditis elegans Natural Diversity Resource](https://www.elegansvariation.org/) for this study. One sequencing run of the reference strain, N2 (SRA accession = DRR142768)⁠ was also used. Due to the large sizes of these files, including them here is not possible.
