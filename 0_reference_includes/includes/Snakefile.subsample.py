SAM_ALIGNERS = ["ngmlr","minimap2"] # Aligners that output sam files
SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]

# Note: Uncomment the following if all isolates will be included in the pipeline
NGMLRDICT10X = {"DL238": "0.06", "N2": "0.13", "ECA36": "0.05", "ECA396": "0.07", "EG4725": "0.05", "JU1400": "0.09", "JU2526": "0.08", "JU2600": "0.06", "JU310": "0.06", "MY2147":"0.05", "MY2693": "0.08", "NIC2": "0.08", "NIC526": "0.07", "QX1794": "0.08","XZ1516": "0.13"}
MINIMAPDICT10X = {"DL238": "0.05", "N2": "0.12", "ECA36": "0.04", "ECA396": "0.07", "EG4725": "0.04", "JU1400": "0.08", "JU2526": "0.08", "JU2600": "0.05", "JU310": "0.05", "MY2147":"0.05", "MY2693": "0.07", "NIC2": "0.08", "NIC526": "0.06", "QX1794": "0.07","XZ1516": "0.12"}
PBMMDICT10X = {"DL238": "0.06", "N2": "0.12", "ECA36": "0.05", "ECA396": "0.07", "EG4725": "0.04", "JU1400": "0.08", "JU2526": "0.08", "JU2600": "0.05", "JU310": "0.06", "MY2147":"0.05", "MY2693": "0.07", "NIC2": "0.08", "NIC526": "0.06", "QX1794": "0.08","XZ1516": "0.13"}
NGMLRDICT20X = {"DL238": "0.12", "N2": "0.26", "ECA36": "0.10", "ECA396": "0.15", "EG4725": "0.09", "JU1400": "0.17", "JU2526": "0.16", "JU2600": "0.11", "JU310": "0.12", "MY2147":"0.11", "MY2693": "0.16", "NIC2": "0.17", "NIC526": "0.13", "QX1794": "0.16","XZ1516": "0.27"}
MINIMAPDICT20X = {"DL238": "0.11", "N2": "0.24", "ECA36": "0.09", "ECA396": "0.14", "EG4725": "0.09", "JU1400": "0.16", "JU2526": "0.16", "JU2600": "0.11", "JU310": "0.11", "MY2147":"0.09", "MY2693": "0.14", "NIC2": "0.15", "NIC526": "0.12", "QX1794": "0.14","XZ1516": "0.23"}
PBMMDICT20X = {"DL238": "0.11", "N2": "0.25", "ECA36": "0.09", "ECA396": "0.15", "EG4725": "0.09", "JU1400": "0.16", "JU2526": "0.16", "JU2600": "0.11", "JU310": "0.12", "MY2147":"0.10", "MY2693": "0.15", "NIC2": "0.16", "NIC526": "0.13", "QX1794": "0.15","XZ1516": "0.26"}
NGMLRDICT40X = {"DL238": "0.24", "N2": "0.52", "ECA36": "0.19", "ECA396": "0.30", "EG4725": "0.18", "JU1400": "0.34", "JU2526": "0.33", "JU2600": "0.22", "JU310": "0.25", "MY2147":"0.21", "MY2693": "0.31", "NIC2": "0.33", "NIC526": "0.27", "QX1794": "0.32","XZ1516": "0.54"}
MINIMAPDICT40X = {"DL238": "0.22", "N2": "0.48", "ECA36": "0.17", "ECA396": "0.27", "EG4725": "0.17", "JU1400": "0.32", "JU2526": "0.31", "JU2600": "0.21", "JU310": "0.22", "MY2147":"0.19", "MY2693": "0.28", "NIC2": "0.30", "NIC526": "0.24", "QX1794": "0.28","XZ1516": "0.47"}
PBMMDICT40X = {"DL238": "0.23", "N2": "0.50", "ECA36": "0.18", "ECA396": "0.29", "EG4725": "0.18", "JU1400": "0.33", "JU2526": "0.32", "JU2600": "0.22", "JU310": "0.23", "MY2147":"0.20", "MY2693": "0.30", "NIC2": "0.31", "NIC526": "0.26", "QX1794": "0.30","XZ1516": "0.51"}
NGMLRDICT60X = {"DL238": "0.36", "N2": "0.78", "ECA36": "0.29", "ECA396": "0.44", "EG4725": "0.28", "JU1400": "0.51", "JU2526": "0.49", "JU2600": "0.34", "JU310": "0.37", "MY2147":"0.32", "MY2693": "0.47", "NIC2": "0.50", "NIC526": "0.40", "QX1794": "0.47","XZ1516": "0.81"}
MINIMAPDICT60X = {"DL238": "0.32", "N2": "0.72", "ECA36": "0.26", "ECA396": "0.41", "EG4725": "0.26", "JU1400": "0.48", "JU2526": "0.47", "JU2600": "0.32", "JU310": "0.33", "MY2147":"0.28", "MY2693": "0.42", "NIC2": "0.45", "NIC526": "0.37", "QX1794": "0.42","XZ1516": "0.70"}
PBMMDICT60X = {"DL238": "0.34", "N2": "0.75", "ECA36": "0.28", "ECA396": "0.44", "EG4725": "0.27", "JU1400": "0.49", "JU2526": "0.49", "JU2600": "0.33", "JU310": "0.35", "MY2147":"0.30", "MY2693": "0.45", "NIC2": "0.47", "NIC526": "0.39", "QX1794": "0.46","XZ1516": "0.77"}

# SUBSAMPLE DATA

## Subsampled FASTQ files to obtain desired depths for NGMLR
rule subsample_ngmlr_10X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/ngmlr/10X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: NGMLRDICT10X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

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

rule subsample_ngmlr_40X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/ngmlr/40X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: NGMLRDICT40X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

rule subsample_ngmlr_60X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/ngmlr/60X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: NGMLRDICT60X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

## Subsampled FASTQ files to obtain desired depths for Minimap2
rule subsample_minimap2_10X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/minimap2/10X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: MINIMAPDICT10X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

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

rule subsample_minimap2_40X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/minimap2/40X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: MINIMAPDICT40X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

rule subsample_minimap2_60X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/minimap2/60X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: MINIMAPDICT60X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

## Subsampled FASTQ files to obtain desired depths for pbmm2
rule subsample_pbmm2_10X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/pbmm2/10X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: PBMMDICT10X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

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

rule subsample_pbmm2_40X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/pbmm2/40X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: PBMMDICT40X[wcs.strain]
	conda:  "../../yaml/seqtk.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"""seqtk sample -s 123 {input} {params.value} > {output}"""

rule subsample_pbmm2_60X:
	input:
		"1_fq_processing/{strain}/original/{strain}_original.fastq"
	output:
		"1_fq_processing/subsampled/pbmm2/60X/{strain}/original/{strain}_original.fastq"
	params:
		value=lambda wcs: PBMMDICT60X[wcs.strain]
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
		"cat {input} | scripts/1_process_input_data/seq-shuf > {output}"

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
