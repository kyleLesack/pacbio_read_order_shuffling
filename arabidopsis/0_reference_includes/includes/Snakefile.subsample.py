SAM_ALIGNERS = ["ngmlr","minimap2"] # Aligners that output sam files
SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]

NGMLRDICT20X = {"1254": "0.50", "6021": "0.26", "6024": "0.36", "9470": "0.31"}
MINIMAPDICT20X = {"1254": "0.56", "6021": "0.28", "6024": "0.39", "9470": "0.34"}
PBMMDICT20X = {"1254": "0.52", "6021": "0.28", "6024": "0.38", "9470": "0.33"}

# SUBSAMPLE BAMS

## Subsampled FASTQ files to obtain 20X depths for NGMLR

rule subsample_ngmlr_20X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/ngmlr/20X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: NGMLRDICT20X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

## Subsampled FASTQ files to obtain 20X depths for Minimap2

rule subsample_minimap2_20X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/minimap2/20X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: MINIMAPDICT20X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

## Subsampled FASTQ files to obtain 20X depths for pbmm2

rule subsample_pbmm2_20X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/pbmm2/20X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: PBMMDICT20X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

# Create FASTQ files with randomized read orders from shuffled original FASTQ files
rule shuffle_subsampled:
	input:
		"1_fq_processing/subsampled/{aligner}/{depth}/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/{aligner}/{depth}/{strain}/{shuffledir}/{strain}_{shuffledir}.fastq"
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
	shell:
		"cat {input} | ../scripts/1_process_input_data/seq-shuf > {output}"

# Align subsampled FASTQ files using pbmm2
rule pbmm2_subsampled:
	input:
		"1_fq_processing/subsampled/pbmm2/{depth}/{strain}/{readorder}/{strain}_{readorder}.fastq"
	output:
		"2_alignments/subsampled/{depth}/pbmm2/{strain}/{readorder}/{strain}_{readorder}_sorted.bam"
	conda:  "../../yaml/pbmm2_1.12.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="06:00:00"
	shell:
		"pbmm2 align  -j {threads} {REFERENCE} {input} {output} --sort --median-filter --sample {wildcards.strain}  --rg '@RG\tID:myid\tSM:{wildcards.strain}'"

# Align subsampled FASTQ files using NGMLR
rule ngmlr_subsampled:
	input:
		"1_fq_processing/subsampled/ngmlr/{depth}/{strain}/{readorder}/{strain}_{readorder}.fastq"
	output:
		temp("2_alignments/subsampled/{depth}/ngmlr/{strain}/{readorder}/{strain}_{readorder}.sam")
	conda:  "../../yaml/ngmlr.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="60:00:00"
	shell:
		"ngmlr -t {threads} -r {REFERENCE} -q {input} -o {output}"

# Align subsampled FASTQ files using Minimap2
rule minimap2_subsampled:
	input:
		"1_fq_processing/subsampled/minimap2/{depth}/{strain}/{readorder}/{strain}_{readorder}.fastq"
	output:
		temp("2_alignments/subsampled/{depth}/minimap2/{strain}/{readorder}/{strain}_{readorder}.sam")
	conda:  "../../yaml/minimap2_2.26.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="4-0:00:00"
	shell:
		"minimap2 -t {threads} -Y -ax map-pb {REFERENCE} {input} > {output}"

# Convert subsampled SAM files to BAM format
rule sam2bam_subsampled:
	input:
		"2_alignments/subsampled/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}.sam"
	output:
		temp("2_alignments/subsampled/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}.bam")
	wildcard_constraints:
		aligner = "|".join(SAM_ALIGNERS)
	conda:  "../../yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 500 + ((attempt - 1) * 10000),
		time_hms="03:00:00"
	shell:
		"samtools view -@ 8 -S -b {input} > {output}"

# Sort subsampled BAM files
rule sortbam_subsampled:
	input:
		"2_alignments/subsampled/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}.bam"
	output:
		"2_alignments/subsampled/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam"
	wildcard_constraints:
		aligner = "|".join(SAM_ALIGNERS)
	conda:  "../../yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="17:00:00"
	shell:
		"samtools sort -o {output} {input} -@ 8"

# Index subsampled BAM files
rule indexbam_subsampled:
	input:
		"2_alignments/subsampled/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam"
	output:
		"2_alignments/subsampled/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam.bai"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS)
	conda:  "../../yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 500 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"samtools index -@ {threads} {input} {output}"
