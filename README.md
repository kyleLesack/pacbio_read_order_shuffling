# Overview

These scripts were used to evaluate the impact of FASTQ file read order on structural variation (SV) calling from long-read DNA sequencing data.

# Requirements

Note: The ./yaml directory include conda yaml environment files. I recommend using conda to install the required software using these files. The required tools and the versions we tested are listed below:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (v.7.25.0)
* [Picard](https://broadinstitute.github.io/picard/) (v.2.27.5)
* [SAMtools](http://www.htslib.org/) (v1.9)
	* Note: SAMTools >= v1.10 didn't work with the NGMLR alignments
* [pybedtools](https://daler.github.io/pybedtools/) (v.0.9.0)
* [seq-shuff](https://github.com/thackl/seq-scripts/blob/master/bin/seq-shuf)
* Aligners
	* [Minimap2](https://github.com/lh3/minimap2) (v.2.26)
	* [NGMLR](https://github.com/philres/ngmlr) (v.0.2.7)
	* [pbmm2](https://github.com/PacificBiosciences/pbmm2) (v.1.12.0)
* SV Callers
	* [pbsv](https://github.com/PacificBiosciences/pbsv) (v2.9.0)
	* [Sniffles](https://github.com/fritzsedlazeck/Sniffles) (v2.2.0)
	* [SVIM](https://github.com/eldariont/svim) (v2.0.0)

## Data Not Included Here

Several files were too large to include here and need to be downloaded or created.

* The [blasr](https://manpages.debian.org/testing/blasr/sawriter.1.en.html) sawriter command was used to generate suffix array files and NGMLR indexes that were required for pbmm2/pbsv
	* The Snakemake pipeline expects them to be in the following directories:
		* 0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa-enc.2.ngm
		* 0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa-ht-13-2.2.ngm
		* 0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa.sa
* *C. elegans PacBio* sequencing data
	* [PacBio sequencing data](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=12908562) from the [Caenorhabditis elegans Natural Diversity Resource](https://www.elegansvariation.org/)
	* [PacBio sequencing data](https://www.ncbi.nlm.nih.gov/sra/?term=DRR142768) for the N2 reference strain

# Snakemake pipeline

The ./Snakefile contains the code to run the pipeline. A second Snakefile contains the instructions required to perform the subsampling steps.

The pipeline expects the input FASTQ files in the following locations:

1_fq_processing/N2/original/N2_original.fastq
1_fq_processing/JU1400/original/JU1400_original.fastq
1_fq_processing/NIC2/original/NIC2_original.fastq
1_fq_processing/JU2526/original/JU2526_original.fastq
1_fq_processing/XZ1516/original/XZ1516_original.fastq
1_fq_processing/MY2693/original/MY2693_original.fastq
1_fq_processing/QX1794/original/QX1794_original.fastq
1_fq_processing/NIC526/original/NIC526_original.fastq
1_fq_processing/DL238/original/DL238_original.fastq
1_fq_processing/ECA396/original/ECA396_original.fastq
1_fq_processing/JU2600/original/JU2600_original.fastq
1_fq_processing/ECA36/original/ECA36_original.fastq
1_fq_processing/EG4725/original/EG4725_original.fastq
1_fq_processing/JU310/original/JU310_original.fastq
1_fq_processing/MY2147/original/MY2147_original.fastq
1_fq_processing/N2/original/N2_original.fastq
