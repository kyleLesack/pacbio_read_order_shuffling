# Overview

These scripts were used to evaluate the impact of FASTQ file read order on structural variation (SV) calling from long-read DNA sequencing data.

# Requirements

Note: The ./yaml directory include conda yaml environment files. I recommend using conda to install the required software using these files. The required tools and the versions we tested are listed below:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (v.7.25.0)
* [Picard Tools](https://broadinstitute.github.io/picard/) (v.2.27.5)
* [SAMtools](http://www.htslib.org/) (v1.9)
	* Note: SAMTools >= v1.10 didn't work with the NGMLR alignments
* [pybedtools](https://daler.github.io/pybedtools/) (v.0.9.0)
* [seq-shuff](https://github.com/thackl/seq-scripts/blob/master/bin/seq-shuf)
* Aligners
	* [Minimap2](https://github.com/lh3/minimap2) (v.2.17)
	* [NGMLR](https://github.com/philres/ngmlr) (v.0.2.7)
	* [pbmm2](https://github.com/PacificBiosciences/pbmm2) (v.1.7.0)
* SV Callers
	* [pbsv](https://github.com/PacificBiosciences/pbsv) (v2.8.0)
	* [Sniffles](https://github.com/fritzsedlazeck/Sniffles) (v2.0.6)
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

The ./Snakefile contains the code to run the pipeline. The pipeline includes several computational steps that were performed on a high performance computing cluster. The default resource requirements that we provided to our SLURM scheduler are specified in the ./slurm/config.yaml file. The default values were overwritten for certain steps, which are specified in certain Snakefile rules under threads and resources.
