include: "0_input/includes/Snakefile.no_shuffle"

ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]
REFERENCE = "/work/wasmuth_lab/mrkyle/sv_calling_pipeline2/1_prepare_reference/output/c_elegans.PRJNA13758.WS263/c_elegans.PRJNA13758.WS263.genomic.fa"
REFERENCE_SAW = "/work/wasmuth_lab/mrkyle/sv_calling_pipeline2/1_prepare_reference/output/c_elegans.PRJNA13758.WS263.sawriter/c_elegans.PRJNA13758.WS263.genomic.fa"
SAM_ALIGNERS = ["ngmlr","minimap2"] # Aligners that output sam files
SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["pbmm2"]
ALL_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]

NGMLRDICT10X = {"DL238": "0.06", "DRR142768": "0.13", "ECA36": "0.05", "ECA396": "0.07", "EG4725": "0.05", "JU1400": "0.09", "JU2526": "0.08", "JU2600": "0.06", "JU310": "0.06", "MY2147":"0.05", "MY2693": "0.08", "NIC2": "0.08", "NIC526": "0.07", "QX1794": "0.08","XZ1516": "0.13"}
MINIMAPDICT10X = {"DL238": "0.05", "DRR142768": "0.12", "ECA36": "0.04", "ECA396": "0.07", "EG4725": "0.04", "JU1400": "0.08", "JU2526": "0.08", "JU2600": "0.05", "JU310": "0.05", "MY2147":"0.05", "MY2693": "0.07", "NIC2": "0.08", "NIC526": "0.06", "QX1794": "0.07","XZ1516": "0.12"}
PBMMDICT10X = {"DL238": "0.06", "DRR142768": "0.12", "ECA36": "0.05", "ECA396": "0.07", "EG4725": "0.04", "JU1400": "0.08", "JU2526": "0.08", "JU2600": "0.05", "JU310": "0.06", "MY2147":"0.05", "MY2693": "0.07", "NIC2": "0.08", "NIC526": "0.06", "QX1794": "0.08","XZ1516": "0.13"}
NGMLRDICT20X = {"DL238": "0.12", "DRR142768": "0.26", "ECA36": "0.10", "ECA396": "0.15", "EG4725": "0.09", "JU1400": "0.17", "JU2526": "0.16", "JU2600": "0.11", "JU310": "0.12", "MY2147":"0.11", "MY2693": "0.16", "NIC2": "0.17", "NIC526": "0.13", "QX1794": "0.16","XZ1516": "0.27"}
MINIMAPDICT20X = {"DL238": "0.11", "DRR142768": "0.24", "ECA36": "0.09", "ECA396": "0.14", "EG4725": "0.09", "JU1400": "0.16", "JU2526": "0.16", "JU2600": "0.11", "JU310": "0.11", "MY2147":"0.09", "MY2693": "0.14", "NIC2": "0.15", "NIC526": "0.12", "QX1794": "0.14","XZ1516": "0.23"}
PBMMDICT20X = {"DL238": "0.11", "DRR142768": "0.25", "ECA36": "0.09", "ECA396": "0.15", "EG4725": "0.09", "JU1400": "0.16", "JU2526": "0.16", "JU2600": "0.11", "JU310": "0.12", "MY2147":"0.10", "MY2693": "0.15", "NIC2": "0.16", "NIC526": "0.13", "QX1794": "0.15","XZ1516": "0.26"}
NGMLRDICT40X = {"DL238": "0.24", "DRR142768": "0.52", "ECA36": "0.19", "ECA396": "0.30", "EG4725": "0.18", "JU1400": "0.34", "JU2526": "0.33", "JU2600": "0.22", "JU310": "0.25", "MY2147":"0.21", "MY2693": "0.31", "NIC2": "0.33", "NIC526": "0.27", "QX1794": "0.32","XZ1516": "0.54"}
MINIMAPDICT40X = {"DL238": "0.22", "DRR142768": "0.48", "ECA36": "0.17", "ECA396": "0.27", "EG4725": "0.17", "JU1400": "0.32", "JU2526": "0.31", "JU2600": "0.21", "JU310": "0.22", "MY2147":"0.19", "MY2693": "0.28", "NIC2": "0.30", "NIC526": "0.24", "QX1794": "0.28","XZ1516": "0.47"}
PBMMDICT40X = {"DL238": "0.23", "DRR142768": "0.50", "ECA36": "0.18", "ECA396": "0.29", "EG4725": "0.18", "JU1400": "0.33", "JU2526": "0.32", "JU2600": "0.22", "JU310": "0.23", "MY2147":"0.20", "MY2693": "0.30", "NIC2": "0.31", "NIC526": "0.26", "QX1794": "0.30","XZ1516": "0.51"}
NGMLRDICT60X = {"DL238": "0.36", "DRR142768": "0.78", "ECA36": "0.29", "ECA396": "0.44", "EG4725": "0.28", "JU1400": "0.51", "JU2526": "0.49", "JU2600": "0.34", "JU310": "0.37", "MY2147":"0.32", "MY2693": "0.47", "NIC2": "0.50", "NIC526": "0.40", "QX1794": "0.47","XZ1516": "0.81"}
MINIMAPDICT60X = {"DL238": "0.32", "DRR142768": "0.72", "ECA36": "0.26", "ECA396": "0.41", "EG4725": "0.26", "JU1400": "0.48", "JU2526": "0.47", "JU2600": "0.32", "JU310": "0.33", "MY2147":"0.28", "MY2693": "0.42", "NIC2": "0.45", "NIC526": "0.37", "QX1794": "0.42","XZ1516": "0.70"}
PBMMDICT60X = {"DL238": "0.34", "DRR142768": "0.75", "ECA36": "0.28", "ECA396": "0.44", "EG4725": "0.27", "JU1400": "0.49", "JU2526": "0.49", "JU2600": "0.33", "JU310": "0.35", "MY2147":"0.30", "MY2693": "0.45", "NIC2": "0.47", "NIC526": "0.39", "QX1794": "0.46","XZ1516": "0.77"}
TEST_STRAINS = ["DL238"]
FULL_SUBSAMPLED_DEPTHS = ["full_depth", "subsampled/10X", "subsampled/20X", "subsampled/40X", "subsampled/60X"]
SUBSAMPLED_DEPTHS = ["10X", "20X", "40X", "60X"]


rule all:
	input:
		#expand("1_fq_processing/shuffled/{strain}/{replicate}/{strain}_shuffled.fastq", strain=TEST_STRAINS, replicate=[1]),
		expand("1_fq_processing/shuffled/{strain}/{replicate}/{strain}_shuffled.fastq", strain=ALL_STRAINS, replicate=[1]),
		expand("2_alignments/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam.bai", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/{aligner}/{strain}/ORIGINAL/{strain}_depth.txt", aligner=SAM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/subsampled/10X/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/subsampled/10X/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("2_alignments/subsampled/20X/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/subsampled/20X/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("2_alignments/subsampled/40X/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/subsampled/40X/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("2_alignments/subsampled/60X/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/subsampled/60X/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam",aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam.bai", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_depth.txt", aligner=SAM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/full_depth/{aligner}/svim/{strain}/{replicate}/variants.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/full_depth/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS),
		expand("3_variant_calls/full_depth/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("3_variant_calls/full_depth/{aligner}/svim/{strain}/ORIGINAL/variants.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{replicate}/variants.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/ORIGINAL/variants.vcf", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60]),
		expand("3_variant_calls/full_depth/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_total_svs.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1],svimqual=[0,15]),
		expand("3_variant_calls/full_depth/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_breakpoints.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1],svimqual=[0,15]),
		expand("3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_total_svs.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1],svimqual=[0,15]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_breakpoints.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1],svimqual=[0,15]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv", aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1]),
		expand("3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, depth=[10,20,40,60], replicate=[1]),
		expand("4_results/{depth}/{analysis}/svim/qual_15/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/{depth}/{analysis}/sniffles/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/{depth}/{analysis}/pbsv/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS),
		expand("4_results/{depth}/{analysis}/svim/qual_15/totals/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/{depth}/{analysis}/sniffles/totals/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/{depth}/{analysis}/pbsv/totals/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS),
		expand("4_results/{depth}/{analysis}/svim/qual_15/ggplot/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/{depth}/{analysis}/sniffles/ggplot/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/{depth}/{analysis}/pbsv/ggplot/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS),
		"5_plots/full_depth/1_sv_agreement.png",
		"5_plots/full_depth/2_breakpoint_agreement-counts.png",
		"5_plots/full_depth/2_breakpoint_agreement.png",
		"5_plots/full_depth/1_sv_agreement-counts.png",
		"5_plots/subsampled/3_sv_agreement_sniffles_svim_subsampled.png",
		"5_plots/subsampled/4_breakpoint_agreement_sniffles_svim_pbsv_subsampled.png",
		"5_plots/subsampled/3_sv_agreement_pbsv_sniffles_svim_subsampled.png",
		"5_plots/subsampled/4_breakpoint_agreement_sniffles_svim_pbsv_subsampled-counts.png",
		"5_plots/subsampled/3_sv_agreement_pbsv_sniffles_svim_subsampled_counts.png"
rule shuffle:
	input:
		"0_input/data/{strain}_all_reads.fastq"
	output:
		"1_fq_processing/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
	shell:
		"cat {input} | scripts/seq-shuf > {output}"

rule pbmm2:
	input:
		"1_fq_processing/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
	output:
		"2_alignments/pbmm2/{strain}/{REP}/{strain}_shuffled.bam"
	conda:  "yaml/pbmm2.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="36:00:00"
	shell:
		"pbmm2 align  -j {threads} {REFERENCE_SAW} {input} {output} --sort --median-filter --sample {wildcards.strain}  --rg '@RG\tID:myid\tSM:EG4725'"

rule ngmlr:
	input:
		"1_fq_processing/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
	output:
		temp("2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled.sam")
	conda:  "yaml/ngmlr.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="36:00:00"
	shell:
		"ngmlr -t {threads} -r {REFERENCE} -q {input} -o {output}"

rule minimap2:
	input:
		"1_fq_processing/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
	output:
		temp("2_alignments/minimap2/{strain}/{REP}/{strain}_shuffled.sam")
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
	wildcard_constraints:
		aligner = "|".join(SAM_ALIGNERS)
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
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="17:00:00"
	shell:
		"samtools sort -o {output} {input} -@ 8"

rule getbamdepth:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_depth.txt"
	wildcard_constraints:
		aligner = "|".join(SAM_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="04:00:00"
	shell:
		"""samtools depth {input} | awk '{{sum+=$3}} END {{ print "Average = ",sum/NR}}' > {output}"""

rule getbamdepthpbmm2:
	input:
		"2_alignments/pbmm2/{strain}/{REP}/{strain}_shuffled.bam"
	output:
		"2_alignments/pbmm2/{strain}/{REP}/{strain}_depth.txt"
	conda:  "yaml/samtools_1.9.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="04:00:00"
	shell:
		"""samtools depth {input} | awk '{{sum+=$3}} END {{ print "Average = ",sum/NR}}' > {output}"""

## SUBSAMPLE BAMS
rule subsamplengmlrto10x:
	input:
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/10X/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: NGMLRDICT10X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsamplengmlrto20x:
	input:
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/20X/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: NGMLRDICT20X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsamplengmlrto40x:
	input:
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/40X/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: NGMLRDICT40X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsamplengmlrto60x:
	input:
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/60X/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: NGMLRDICT60X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsampleminimap2to10x:
	input:
		"2_alignments/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/10X/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: MINIMAPDICT10X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsampleminimap2to20x:
	input:
		"2_alignments/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/20X/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: MINIMAPDICT20X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsampleminimap2to40x:
	input:
		"2_alignments/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/40X/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: MINIMAPDICT40X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsampleminimap2to60x:
	input:
		"2_alignments/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/60X/minimap2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: MINIMAPDICT60X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsamplepbmm2to10x:
	input:
		"2_alignments/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/10X/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: PBMMDICT10X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsamplepbmm2to20x:
	input:
		"2_alignments/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/20X/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: PBMMDICT20X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsamplepbmm2to40x:
	input:
		"2_alignments/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/40X/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: PBMMDICT40X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule subsamplepbmm2to60x:
	input:
		"2_alignments/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/60X/pbmm2/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: PBMMDICT60X[wcs.strain]
	wildcard_constraints:
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

rule pbsvdiscover:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled.bam"
	output:
		temp("3_variant_calls/full_depth/{aligner}/{strain}/{REP}/{strain}.svsig.gz")
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="02:00:00"
	shell:
		"pbsv discover {input} {output}"

rule pbsvcall:
	input:
		"3_variant_calls/full_depth/{aligner}/{strain}/{REP}/{strain}.svsig.gz"
	output:
		"3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{REP}/{strain}.vcf"
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="08:00:00"
	shell:
		"pbsv call -j 8 {REFERENCE_SAW} {input} {output}"

rule svim:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"3_variant_calls/full_depth/{aligner}/svim/{strain}/{REP}/variants.vcf"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	params:
		outdir="3_variant_calls/full_depth/{aligner}/svim/{strain}/{REP}/"
	conda:  "yaml/svim.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"svim alignment {params.outdir} {input} {REFERENCE}"

rule indexbam:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam.bai"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS)
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"samtools index -@ {threads} {input} {output}"

rule sniffles:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam",
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam.bai"
	output:
		"3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{REP}/{strain}.vcf"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	params:
		bamfile="2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	conda:  "yaml/sniffles.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"sniffles --input {params.bamfile} -t 8 --vcf {output}"

rule compare_svim:
	input:
		"3_variant_calls/full_depth/{aligner}/svim/{strain}/{replicate}/variants.vcf",
		"3_variant_calls/full_depth/{aligner}/svim/{strain}/ORIGINAL/variants.vcf"
	output:
		"3_variant_calls/full_depth/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_total_svs.csv",
		"3_variant_calls/full_depth/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
		svimqual = "|".join(["0","15"])
	params:
		original = "3_variant_calls/full_depth/{aligner}/svim/{strain}/ORIGINAL/variants.vcf",
		shuffled = "3_variant_calls/full_depth/{aligner}/svim/{strain}/{replicate}/variants.vcf",
		quality = "{svimqual}"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.original} svim --minsize 0 --min_qual_svim {params.quality}"

rule compare_sniffles:
	input:
		"3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf",
		"3_variant_calls/full_depth/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf"
	output:
		"3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv",
		"3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
	params:
		original = "3_variant_calls/full_depth/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf",
		shuffled = "3_variant_calls/full_depth/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.original} sniffles --minsize 0"

rule compare_pbsv:
	input:
		"3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf",
		"3_variant_calls/full_depth/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf"
	output:
		"3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv",
		"3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
	params:
		original = "3_variant_calls/full_depth/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf",
		shuffled = "3_variant_calls/full_depth/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.original} pbsv --minsize 0"

rule concatenate_total_svs_sniffles:
	input:
		expand("3_variant_calls/{{aligner}}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/{aligner}/sniffles/total_svs.csv"
	params:
		cmd="tail -q -n 1",
		outdir = "4_summaries/concatenated/{aligner}/sniffles/"
	shell:
		"mkdir -p {params.outdir} && {params.cmd} {input} > {output} "

rule concatenate_total_breakpoints_sniffles:
	input:
		expand("3_variant_calls/{{aligner}}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/{aligner}/sniffles/total_breakpoints.csv"
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/concatenated/{aligner}/sniffles/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_svs_svim:
	input:
		expand("3_variant_calls/{{aligner}}/svim/{strain}/{replicate}/QUAL_0/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/{aligner}/svim/QUAL_0/total_svs.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/concatenated/{aligner}/svim/QUAL_0/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_breakpoints_svim:
	input:
		expand("3_variant_calls/{{aligner}}/svim/{strain}/{replicate}/QUAL_0/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/{aligner}/svim/QUAL_0/total_breakpoints.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/concatenated/{aligner}/svim/QUAL_0/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_svs_pbsv:
	input:
		expand("3_variant_calls/pbmm2/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/pbmm2/pbsv/total_svs.csv",
	params:
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/concatenated/pbmm2/pbsv/"
	shell:
		"scripts/concatenate_csvs.sh '{input}' {params.outdir} {output}"

rule concatenate_total_breakpoints_pbsv:
	input:
		expand("3_variant_calls/pbmm2/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/pbmm2/pbsv/total_breakpoints.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/concatenated/pbmm2/pbsv/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule pbsvdiscover_subsampled:
	input:
		"2_alignments/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		temp("3_variant_calls/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}.svsig.gz")
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"]),
		depth = "|".join(["10","20","40","60"])
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="02:00:00"
	shell:
		"pbsv discover {input} {output}"

rule pbsvcall_subsampled:
	input:
		"3_variant_calls/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}.svsig.gz"
	output:
		"3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{REP}/{strain}.vcf"
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"]),
		depth = "|".join(["10","20","40","60"])
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="08:00:00"
	shell:
		"pbsv call -j 8 {REFERENCE_SAW} {input} {output}"

rule svim_subsampled:
	input:
		"2_alignments/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{REP}/variants.vcf"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"]),
		depth = "|".join(["10","20","40","60"])
	params:
		outdir="3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{REP}/"
	conda:  "yaml/svim.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"svim alignment {params.outdir} {input} {REFERENCE}"

rule indexbam_subsampled:
	input:
		"2_alignments/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam.bai"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS)
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"samtools index -@ {threads} {input} {output}"

rule sniffles_subsampled:
	input:
		"2_alignments/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam",
		"2_alignments/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam.bai"
	output:
		"3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{REP}/{strain}.vcf"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	params:
		bamfile="2_alignments/subsampled/{depth}X/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	conda:  "yaml/sniffles.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"sniffles --input {params.bamfile} -t 8 --vcf {output}"

rule concatenate_total_svs_sniffles_subsampled:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/{{aligner}}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/{aligner}/sniffles/total_svs.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/{aligner}/sniffles/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_breakpoints_sniffles_subsampled:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/{{aligner}}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/{aligner}/sniffles/total_breakpoints.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/{aligner}/sniffles/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_svs_svim_subsampled:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/{{aligner}}/svim/{strain}/{replicate}/QUAL_0/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_0/total_svs.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_0/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_breakpoints_svim_subsampled:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/{{aligner}}/svim/{strain}/{replicate}/QUAL_0/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_0/total_breakpoints.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_0/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_svs_pbsv_subsampled:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/pbmm2/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/pbmm2/pbsv/total_svs.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/pbmm2/pbsv/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule concatenate_total_breakpoints_pbsv_subsampled:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/pbmm2/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/pbmm2/pbsv/total_breakpoints.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/pbmm2/pbsv/"
	shell:
		r"""mkdir -p {params.outdir} && {params.cmd} {input} > {output}"""

rule compare_svim_subsampled:
	input:
		"3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{replicate}/variants.vcf",
		"3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/ORIGINAL/variants.vcf"
	output:
		"3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_total_svs.csv",
		"3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{replicate}/QUAL_{svimqual}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
		svimqual = "|".join(["0","15"])
	params:
		original = "3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/ORIGINAL/variants.vcf",
		shuffled = "3_variant_calls/subsampled/{depth}X/{aligner}/svim/{strain}/{replicate}/variants.vcf",
		quality = "{svimqual}"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.original} svim --minsize 0 --min_qual_svim {params.quality}"

rule compare_sniffles_subsampled:
	input:
		"3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf",
		"3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf"
	output:
		"3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv",
		"3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
	params:
		original = "3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf",
		shuffled = "3_variant_calls/subsampled/{depth}X/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.original} sniffles --minsize 0"

rule compare_pbsv_subsampled:
	input:
		"3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf",
		"3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf"
	output:
		"3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv",
		"3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
	params:
		original = "3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf",
		shuffled = "3_variant_calls/subsampled/{depth}X/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.original} pbsv --minsize 0"

rule concatenate_total_svs_svim_qual15:
	input:
		expand("3_variant_calls/{{aligner}}/svim/{strain}/{replicate}/QUAL_15/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/{aligner}/svim/QUAL_15/total_svs.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/concatenated/{aligner}/svim/QUAL_15/"
	shell:
		r"""scripts/concatenate_csvs.sh "{input}" {params.outdir} > {output}"""

rule concatenate_total_breakpoints_svim_qual15:
	input:
		expand("3_variant_calls/{{aligner}}/svim/{strain}/{replicate}/QUAL_15/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/concatenated/{aligner}/svim/QUAL_15/total_breakpoints.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/concatenated/{aligner}/svim/QUAL_15/"
	shell:
		r"""scripts/concatenate_csvs.sh "{input}" {params.outdir} > {output}"""

rule concatenate_total_svs_svim_subsampled_qual15:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/{{aligner}}/svim/{strain}/{replicate}/QUAL_15/summary/summary_total_svs.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_15/total_svs.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_15/"
	shell:
		r"""scripts/concatenate_csvs.sh {input} {params.outdir} > {output} """

rule concatenate_total_breakpoints_svim_subsampled_qual15:
	input:
		expand("3_variant_calls/subsampled/{{depth}}X/{{aligner}}/svim/{strain}/{replicate}/QUAL_15/summary/summary_breakpoints.csv", strain=ALL_STRAINS, replicate=[1]),
	output:
		"4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_15/total_breakpoints.csv",
	params:
		cmd="tail -q -n 1",
		outdir="/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/4_summaries/subsampled/{depth}X/concatenated/{aligner}/svim/QUAL_15/"
	shell:
		r"""scripts/concatenate_csvs.sh {input} {params.outdir} > {output}"""

# Summarize results by type
rule summarize_svim_results_by_type:
	input:
		expand("3_variant_calls/{depth}/{aligner}/svim/{strain}/1/QUAL_15/summary/summary_total_svs.csv", depth = FULL_SUBSAMPLED_DEPTHS, aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS),
		expand("3_variant_calls/{aligner}/svim/{strain}/1/QUAL_15/summary/summary_breakpoints.csv", depth = FULL_SUBSAMPLED_DEPTHS, aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS),
	output:
		expand("4_results/{depth}/{analysis}/svim/qual_15/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/1_summarize_by_sv_type.py"

rule summarize_sniffles_results_by_type:
	input:
		expand("3_variant_calls/{depth}/{aligner}/sniffles/{strain}/1/summary/summary_total_svs.csv", depth = FULL_SUBSAMPLED_DEPTHS, aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS),
		expand("3_variant_calls/{aligner}/sniffles/{strain}/1/summary/summary_breakpoints.csv", depth = FULL_SUBSAMPLED_DEPTHS, aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS),
	output:
		expand("4_results/{depth}/{analysis}/sniffles/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/1_summarize_by_sv_type.py"

rule summarize_pbsv_results_by_type:
	input:
		expand("3_variant_calls/{depth}/{aligner}/pbsv/{strain}/1/summary/summary_total_svs.csv", depth = FULL_SUBSAMPLED_DEPTHS, aligner = PBSV_ALIGNERS, strain = ALL_STRAINS),
		expand("3_variant_calls/{aligner}/pbsv/{strain}/1/summary/summary_breakpoints.csv", depth = FULL_SUBSAMPLED_DEPTHS, aligner = PBSV_ALIGNERS, strain = ALL_STRAINS),
	output:
		expand("4_results/{depth}/{analysis}/pbsv/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS, strain = ALL_STRAINS)
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/1_summarize_by_sv_type.py"

rule summarize_svim_results_by_type_add_totals:
	input:
		expand("4_results/{depth}/{analysis}/svim/qual_15/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	output:
		expand("4_results/{depth}/{analysis}/svim/qual_15/totals/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/2_add_totals.py"

rule summarize_sniffles_results_by_type_add_totals:
	input:
		expand("4_results/{depth}/{analysis}/sniffles/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	output:
		expand("4_results/{depth}/{analysis}/sniffles/totals/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/2_add_totals.py"

rule summarize_pbsv_results_by_type_add_totals:
	input:
		expand("4_results/{depth}/{analysis}/pbsv/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS, strain = ALL_STRAINS)
	output:
		expand("4_results/{depth}/{analysis}/pbsv/totals/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS, strain = ALL_STRAINS)
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/2_add_totals.py"

rule summarize_svim_results_by_type_create_data4ggplot:
	input:
		expand("4_results/{depth}/{analysis}/svim/qual_15/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	output:
		expand("4_results/{depth}/{analysis}/svim/qual_15/ggplot/svim-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/3_create_ggplot_data.py"

rule summarize_sniffles_results_by_type_create_data4ggplot:
	input:
		expand("4_results/{depth}/{analysis}/sniffles/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	output:
		expand("4_results/{depth}/{analysis}/sniffles/ggplot/sniffles-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS)
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/3_create_ggplot_data.py"

rule summarize_pbsv_results_by_type_create_data4ggplot:
	input:
		expand("4_results/{depth}/{analysis}/pbsv/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS, strain = ALL_STRAINS)
	output:
		expand("4_results/{depth}/{analysis}/pbsv/ggplot/pbsv-{aligner}.csv", depth = FULL_SUBSAMPLED_DEPTHS, analysis = ["sv_intersection_agreement", "breakpoint_agreement"], aligner = PBSV_ALIGNERS, strain = ALL_STRAINS)
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/summarize_results_by_type/3_create_ggplot_data.py"

rule r_plot_sv_agreement:
	input:
		"4_results/full_depth/sv_intersection_agreement/pbsv/ggplot/pbsv-pbmm2.csv",
		"4_results/full_depth/sv_intersection_agreement/sniffles/ggplot/sniffles-minimap2.csv",
		"4_results/full_depth/sv_intersection_agreement/sniffles/ggplot/sniffles-ngmlr.csv",
		"4_results/full_depth/sv_intersection_agreement/sniffles/ggplot/sniffles-pbmm2.csv",
		"4_results/full_depth/sv_intersection_agreement/svim/qual_15/ggplot/svim-minimap2.csv",
		"4_results/full_depth/sv_intersection_agreement/svim/qual_15/ggplot/svim-ngmlr.csv",
		"4_results/full_depth/sv_intersection_agreement/svim/qual_15/ggplot/svim-pbmm2.csv"	
	output:
		"5_plots/full_depth/1_sv_agreement.png",
		"5_plots/full_depth/1_sv_agreement-counts.png"
	conda:  "yaml/R.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 1000),
		time="00:30:00"
	script:
		"scripts/r/1_sv_agreement.R"

rule r_plot_breakpoint_agreement:
	input:
		"4_results/full_depth/breakpoint_agreement/pbsv/ggplot/pbsv-pbmm2.csv",
		"4_results/full_depth/breakpoint_agreement/sniffles/ggplot/sniffles-minimap2.csv",
		"4_results/full_depth/breakpoint_agreement/sniffles/ggplot/sniffles-ngmlr.csv",
		"4_results/full_depth/breakpoint_agreement/sniffles/ggplot/sniffles-pbmm2.csv",
		"4_results/full_depth/breakpoint_agreement/svim/qual_15/ggplot/svim-minimap2.csv",
		"4_results/full_depth/breakpoint_agreement/svim/qual_15/ggplot/svim-ngmlr.csv",
		"4_results/full_depth/breakpoint_agreement/svim/qual_15/ggplot/svim-pbmm2.csv"	
	output:
		"5_plots/full_depth/2_breakpoint_agreement.png",
		"5_plots/full_depth/2_breakpoint_agreement-counts.png"
	conda:  "yaml/R.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 1000),
		time="00:30:00"
	script:
		"scripts/r/2_breakpoint_agreement.R"


rule r_plot_sv_agreement_sniffles_svim:
	input:
		expand("4_results/subsampled/{depth}/sv_intersection_agreement/pbsv/ggplot/pbsv-pbmm2.csv", depth = ["10X", "20X", "40X", "60X"]),
		expand("4_results/subsampled/{depth}/sv_intersection_agreement/sniffles/ggplot/sniffles-{aligner}.csv", depth = ["10X", "20X", "40X", "60X"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/subsampled/{depth}/sv_intersection_agreement/svim/qual_15/ggplot/svim-{aligner}.csv", depth = ["10X", "20X", "40X", "60X"], aligner = SNIFFLES_SVIM_ALIGNERS)
	output:
		"5_plots/subsampled/3_sv_agreement_sniffles_svim_subsampled.png",
		"5_plots/subsampled/3_sv_agreement_pbsv_sniffles_svim_subsampled.png",
		"5_plots/subsampled/3_sv_agreement_pbsv_sniffles_svim_subsampled_counts.png"
	conda:  "yaml/R.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 1000),
		time="00:30:00"
	script:
		"scripts/r/3_sv_agreement_subsampled_pbsv_sniffles_svim.R"

rule r_plot_breakpoint_agreement_pbsv_sniffles_svim_subsampled:
	input:
		expand("4_results/subsampled/{depth}/breakpoint_agreement/pbsv/ggplot/pbsv-pbmm2.csv", depth = ["10X", "20X", "40X", "60X"]),
		expand("4_results/subsampled/{depth}/breakpoint_agreement/sniffles/ggplot/sniffles-{aligner}.csv", depth = ["10X", "20X", "40X", "60X"], aligner = SNIFFLES_SVIM_ALIGNERS),
		expand("4_results/subsampled/{depth}/breakpoint_agreement/svim/qual_15/ggplot/svim-{aligner}.csv", depth = ["10X", "20X", "40X", "60X"], aligner = SNIFFLES_SVIM_ALIGNERS)
	output:
		"5_plots/subsampled/4_breakpoint_agreement_sniffles_svim_pbsv_subsampled.png",
		"5_plots/subsampled/4_breakpoint_agreement_sniffles_svim_pbsv_subsampled-counts.png"
	conda:  "yaml/R.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 1000),
		time="00:30:00"
	script:
		"scripts/r/4_breakpoint_agreement_subsampled_sniffles_svim_pbsv.R"
