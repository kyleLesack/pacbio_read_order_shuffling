SPLIT_STRAINS = ["DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"] # Large strains that took too long to shuffle. I split the Fastq files before shuffling and concatenation.
NO_SPLIT_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768"]
ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]
REFERENCE = "/work/wasmuth_lab/mrkyle/sv_calling_pipeline2/1_prepare_reference/output/c_elegans.PRJNA13758.WS263/c_elegans.PRJNA13758.WS263.genomic.fa"
ALIGNERS = ["ngmlr","minimap2"]
SPLIT_DIRS = ["no_split", "needs_split"]

rule all:
	input:
		expand("1_fq_processing/all_strains/shuffled/{no_split_strain}/{replicate}/{no_split_strain}_shuffled.fastq", no_split_strain=NO_SPLIT_STRAINS, replicate=[1,2,3,4,5]),
		expand("1_fq_processing/all_strains/shuffled/{split_strain}/{replicate}/{split_strain}_shuffled.fastq", split_strain=SPLIT_STRAINS, replicate=[1,2,3,4,5]),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_shuffled.sam", aligner=ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/sniffles/{strain}/{replicate}/{strain}.vcf", strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/svim/{strain}/{replicate}/variants.vcf", strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
rule split:
	input:
		"0_input/data/{split_strain}_all_reads.fastq"
	output:
		"1_fq_processing/needs_split/1_splits/{split_strain}/{split_strain}_all_reads.part_001.fastq",
		"1_fq_processing/needs_split/1_splits/{split_strain}/{split_strain}_all_reads.part_002.fastq"
	conda:  "yaml/seqkit.yaml"
	params:
		outdir="1_fq_processing/needs_split/splits/{split_strain}/"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 1000),
		time_hms="00:30:00"
	shell:
		"seqkit split2 {input} -p 2 -O {params.outdir}"

rule shuffle_split_1:
	input:
		"1_fq_processing/needs_split/1_splits/{split_strain}/{split_strain}_all_reads.part_001.fastq"
	output:
		#temp("1_fq_processing/needs_split/shuffled_split/{strain}/{REP}/{strain}_shuffled.part_001.fastq")
		temp("1_fq_processing/needs_split/2_shuffled_split/{split_strain}/{REP}/{split_strain}_shuffled.part_001.fastq")
	conda:  "yaml/bbtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
	shell:
		"shuffle.sh {input} out={output}"

rule shuffle_split_2:
	input:
		"1_fq_processing/needs_split/1_splits/{split_strain}/{split_strain}_all_reads.part_002.fastq"
	output:
		#temp("1_fq_processing/needs_split/shuffled_split/{strain}/{REP}/{strain}_shuffled.part_002.fastq")
		temp("1_fq_processing/needs_split/2_shuffled_split/{split_strain}/{REP}/{split_strain}_shuffled.part_002.fastq")
	conda:  "yaml/bbtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
	shell:
		"shuffle.sh {input} out={output}"

rule combine_split:
	input:
		split_1="1_fq_processing/needs_split/2_shuffled_split/{split_strain}/{REP}/{split_strain}_shuffled.part_001.fastq",
		split_2="1_fq_processing/needs_split/2_shuffled_split/{split_strain}/{REP}/{split_strain}_shuffled.part_002.fastq"
	output:
		#"1_fq_processing/needs_split/shuffled/{strain}/{REP}/{strain}_shuffled_combined.fastq"
		"1_fq_processing/all_strains/shuffled/{split_strain}/{REP}/{strain}_shuffled.fastq"
	wildcard_constraints:
		split_strain = "|".join(SPLIT_STRAINS)
	resources:
		mem_mb=lambda _, attempt: 2000 + ((attempt - 1) * 5000),
                time_hms="05:00:00"
	shell:
		"cat {input.split_1} {input.split_2} > {output}"

rule shuffle:
	input:
		"0_input/data/{no_split_strain}_all_reads.fastq"
	output:
		#"1_fq_processing/no_split/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
		"1_fq_processing/all_strains/shuffled/{no_split_strain}/{REP}/{no_split_strain}_shuffled.fastq"
	wildcard_constraints:
		no_split_strain = "|".join(NO_SPLIT_STRAINS)        
	conda:  "yaml/bbtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 90000 + ((attempt - 1) * 10000),
	shell:
		"shuffle.sh {input} out={output}"

rule ngmlr:
	input:
		"1_fq_processing/all_strains/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
	output: 
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled.sam"
	conda:  "yaml/ngmlr.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="36:00:00"
	shell:
		"ngmlr -t {threads} -r {REFERENCE} -q {input} -o {output}"

rule minimap2:
	input:
		"1_fq_processing/all_strains/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
	output: 
		"2_alignments/minimap2/{strain}/{REP}/{strain}_shuffled.sam"
	conda:  "yaml/pbhoney.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 100000 + ((attempt - 1) * 10000),
		time="36:00:00"
	shell:
		"minimap2 -t {threads} -ax map-pb {REFERENCE} {input} > {output}"

rule sam2bam:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled.sam"
	output: 
		temp("2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled.bam")
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"samtools view -@ 8 -S -b {input} > {output}"

rule sortbam:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled.bam"
	output:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="17:00:00"
	shell:
		"samtools sort -o {output} {input} -@ 8"

rule svim:
	input:
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"3_variant_calls/svim/{strain}/{REP}/variants.vcf"
	params:
		outdir="3_variant_calls/svim/{strain}/{REP}/"
	conda:  "yaml/svim.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"svim alignment {params.outdir} {input} {REFERENCE}"

rule sniffles:
	input:
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"3_variant_calls/sniffles/{strain}/{REP}/{strain}.vcf"
	conda:  "yaml/sniffles.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"sniffles -m {input} -t 8 -v {output}"

