include: "0_reference_includes/includes/Snakefile.subsample.py"
ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "N2", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]

SHUFFLED_DIRS = ["shuffled1","shuffled2","shuffled3","shuffled4","shuffled5"]
READORDER = ["shuffled1","shuffled2","shuffled3","shuffled4","shuffled5","original"]

REFERENCE = "0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa" # C. elegans reference genome
SAM_ALIGNERS = ["ngmlr","minimap2"] # Aligners that output sam files
SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
SVIM_QUAL =["0"] # SVIM quality filter. Add higher values to list to evaluate how FASTQ read order affects SV calling in SVIM with more stringent filtering threshold(s)

FULL_SUBSAMPLED_DEPTHS = ["full_depth", "subsampled/10X", "subsampled/20X", "subsampled/40X", "subsampled/60X"]
SUBSAMPLED_DEPTHS = ["10X", "20X", "40X", "60X"]

COMPARISON_SUMMARY_FILES = ["summary/overlap_comparison_all_svs.csv", "unique/non_intersecting_all_svs.csv", "summary/overlap_comparison_coords_all_svs.csv", "unique/non_intersecting_coords_all_svs.csv", "summary/overlap_comparison_relaxed_all_svs.csv", "unique/non_intersecting_relaxed_all_svs.csv"]
OVERLAP_SUMMARY_FILES = ['agreement_summary_total.csv', 'agreement_summary_coords_total.csv', 'agreement_summary_relaxed_total.csv']
SUBSAMPLED_COMBINED_TABLE_FILES = ["total.csv", "coords_total.csv", "relaxed_total.csv"]
PLOTS = ["full_depth/1_sv_agreement_grouped_by_caller_counts_top_labels.png", "full_depth/coords/1_sv_agreement_coords_grouped_by_caller_counts_top_labels.png", "subsampled/2_sv_agreement_pbsv_counts.png", "subsampled/2_sv_agreement_sniffles_counts_subsampled.png", "subsampled/2_sv_agreement_pbsv_sniffles_svim_subsampled_counts.png", "subsampled/2_sv_agreement_svim_counts_subsampled.png"]
REPLICATES = ["1","2","3","4","5"]
COMPARISON_FILES = ["results.csv", "results_no_meta.csv", "results_proportions.csv", "results_proportions_no_meta.csv", "sv_intersection_summary.log", "sv_intersection_summary_no_meta.log", "vcf_summary.log"]
SV_CALLERS = ["pbsv", "sniffles", "svim"]
AGREEMENT_SUMMARY_FILES = ["combined_results.csv", "vcf_agreement_results.csv"]
AGREEMENT_SUMMARY_NO_META_FILES = ["combined_results_no_meta.csv", "vcf_agreement_results_no_meta.csv"]
ORIGINAL_COMPARISON_FILES = ["vcf_summary.log", "sv_intersection_stats.log", "sv_intersection_stats_no_meta.log", "results.csv", "results_proportions.csv", "results_no_meta.csv", "results_proportions_no_meta.csv"]

rule all:
	input:
		expand("1_fq_processing/{strain}/{shuffledir}/{strain}_shuffled.fastq", strain=ALL_STRAINS, shuffledir = SHUFFLED_DIRS),
		expand("1_fq_processing/subsampled/{aligner}/{depth}/{strain}/original/{strain}_original.fastq", aligner=SNIFFLES_SVIM_ALIGNERS, depth=SUBSAMPLED_DEPTHS, strain=ALL_STRAINS),
		expand("1_fq_processing/subsampled/{aligner}/{depth}/{strain}/{shuffledir}/{strain}_{shuffledir}.fastq", aligner=SNIFFLES_SVIM_ALIGNERS, depth=SUBSAMPLED_DEPTHS,strain=ALL_STRAINS, shuffledir = SHUFFLED_DIRS),
		#expand("2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.{extension}", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, readorder=READORDER, extension=["bam","bam.bai"]),
		#expand("3_variant_calls/{depth}/{aligner}/pbsv/{strain}/{readorder}/{strain}.vcf", depth=FULL_SUBSAMPLED_DEPTHS, aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, readorder=READORDER),
		#expand("3_variant_calls/{depth}/{aligner}/sniffles/{strain}/{readorder}/{strain}.vcf", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, readorder=READORDER),
		#expand("3_variant_calls/{depth}/{aligner}/svim/{strain}/{readorder}/variants.vcf", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, readorder=READORDER),
		#expand("4_variant_call_comparisons/{depth}/{aligner}/pbsv/{strain}/{results_file}", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, results_file=COMPARISON_FILES),
		#expand("4_variant_call_comparisons/{depth}/{aligner}/sniffles/{strain}/{results_file}", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, results_file=COMPARISON_FILES),
		#expand("4_variant_call_comparisons/{depth}/{aligner}/svim/{strain}/{results_file}", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, strain=ALL_STRAINS, results_file=COMPARISON_FILES),
		#expand("4_variant_call_comparisons/{depth}/{aligner}/{caller}/{strain}/original_shuffled_comparisons/{readorder}/{results_file}", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, caller=SV_CALLERS, strain=ALL_STRAINS, readorder=SHUFFLED_DIRS, results_file=ORIGINAL_COMPARISON_FILES),
		#expand("5_vcf_comparison_summaries/{depth}/{aligner}/{caller}/{results_file}", depth=FULL_SUBSAMPLED_DEPTHS, aligner=SNIFFLES_SVIM_ALIGNERS, caller = SV_CALLERS, results_file = ["csv_totals_all_strains_replicates_jaccard.csv", "csv_totals_all_strains_replicates_differences.csv"]),
		#expand("5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/{strain}_{results_file}", depth = FULL_SUBSAMPLED_DEPTHS, aligner = SNIFFLES_SVIM_ALIGNERS, strain = ALL_STRAINS, results_file = ["results.csv","results_no_meta.csv"]),
		#expand("5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/{results_file}", depth = FULL_SUBSAMPLED_DEPTHS, aligner = SNIFFLES_SVIM_ALIGNERS, results_file = ["results_combined.csv","results_combined_no_meta.csv"]),

# Create FASTQ files with randomized read orders from original FASTQ files
rule shuffle:
	input:
		ancient("1_fq_processing/{strain}/original/{strain}_original.fastq")
	output:
		"1_fq_processing/{strain}/{shuffledir}/{strain}_shuffled.fastq"
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
	shell:
		"cat {input} | scripts/1_process_input_data/seq-shuf > {output}"

# Align FASTQ files using pbmm2
rule pbmm2:
	input:
		ancient("1_fq_processing/{strain}/{readorder}/{strain}_{readorder}.fastq")
	output:
		"2_alignments/full_depth/pbmm2/{strain}/{readorder}/{strain}_{readorder}_sorted.bam"
	conda:  "yaml/pbmm2_1.12.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="06:00:00"
	shell:
		"pbmm2 align  -j {threads} {REFERENCE} {input} {output} --sort --median-filter --sample {wildcards.strain}  --rg '@RG\tID:myid\tSM:{wildcards.strain}'"

# Align FASTQ files using NGMLR
rule ngmlr:
	input:
		ancient("1_fq_processing/{strain}/{readorder}/{strain}_{readorder}.fastq")
	output:
		temp("2_alignments/full_depth/ngmlr/{strain}/{readorder}/{strain}_{readorder}.sam")
	conda:  "yaml/ngmlr.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="60:00:00"
	shell:
		"ngmlr -t {threads} -r {REFERENCE} -q {input} -o {output}"

# Align FASTQ files using Minimap2
rule minimap2:
	input:
		ancient("1_fq_processing/{strain}/{readorder}/{strain}_{readorder}.fastq")
	output:
		temp("2_alignments/full_depth/minimap2/{strain}/{readorder}/{strain}_{readorder}.sam")
	conda:  "yaml/minimap2_2.26.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="4-0:00:00"
	shell:
		"minimap2 -t {threads} -Y -ax map-pb {REFERENCE} {input} > {output}"

# Convert SAM files to BAM format
rule sam2bam:
	input:
		"2_alignments/full_depth/{aligner}/{strain}/{readorder}/{strain}_{readorder}.sam"
	output:
		temp("2_alignments/full_depth/{aligner}/{strain}/{readorder}/{strain}_{readorder}.bam")
	wildcard_constraints:
		aligner = "|".join(SAM_ALIGNERS)
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 500 + ((attempt - 1) * 10000),
		time_hms="03:00:00"
	shell:
		"samtools view -@ 8 -S -b {input} > {output}"

# Sort BAM files
rule sortbam:
	input:
		"2_alignments/full_depth/{aligner}/{strain}/{readorder}/{strain}_{readorder}.bam"
	output:
		"2_alignments/full_depth/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam"
	wildcard_constraints:
		aligner = "|".join(SAM_ALIGNERS)
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="17:00:00"
	shell:
		"samtools sort -o {output} {input} -@ 8"

# Index BAM files
rule indexbam:
	input:
		"2_alignments/full_depth/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam"
	output:
		"2_alignments/full_depth/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam.bai"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS)
	conda:  "yaml/samtools_1.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 500 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
		"samtools index -@ {threads} {input} {output}"

# Discover SV signatures with pbsv
rule pbsvdiscover:
	input:
		bamfile="2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam",
		bamindex="2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam.bai"
	output:
		temp("3_variant_calls/{depth}/{aligner}/pbsv/{strain}/{readorder}/{strain}.svsig.gz")
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS)
	params:
		samplename="sample_{strain}"
	conda:  "yaml/pbsv_2.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 25000 + ((attempt - 1) * 25000),
		time="02:00:00"
	shell:
		"pbsv discover -s {params.samplename} {input.bamfile} {output}"

# Call SVs with pbsv using pbsv signature files
rule pbsvcall:
	input:
		"3_variant_calls/{depth}/{aligner}/pbsv/{strain}/{readorder}/{strain}.svsig.gz"
	output:
		"3_variant_calls/{depth}/{aligner}/pbsv/{strain}/{readorder}/{strain}.vcf"
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS)
	conda:  "yaml/pbsv_2.9.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="08:00:00"
	shell:
		"pbsv call -j 8 {REFERENCE} {input} {output}"

# Call SVs with SVIM
rule svim:
	input:
		bamfile="2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam",
		bamindex="2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam.bai"
	output:
		"3_variant_calls/{depth}/{aligner}/svim/{strain}/{readorder}/variants.vcf"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
	params:
		outdir="3_variant_calls/{depth}/{aligner}/svim/{strain}/{readorder}/"
	conda:  "yaml/svim2.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 20000 + ((attempt - 1) * 10000),
		time_hms="02:00:00"
	shell:
		"svim alignment {params.outdir} {input.bamfile} {REFERENCE}"

# Call SVs with Sniffles
rule sniffles:
	input:
		"2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam",
		"2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam.bai"
	output:
		"3_variant_calls/{depth}/{aligner}/sniffles/{strain}/{readorder}/{strain}.vcf"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_SVIM_ALIGNERS),
	params:
		bamfile="2_alignments/{depth}/{aligner}/{strain}/{readorder}/{strain}_{readorder}_sorted.bam"
	conda:  "yaml/sniffles_2.2.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="02:00:00"
	shell:
		"sniffles --input {params.bamfile} -t 8 --vcf {output}"

# Compare the pbsv SVs predicted from the original FASTQ files with the shuffled versions
rule compare_pbsv:
	input:
		expand("3_variant_calls/{depth}/{aligner}/pbsv/{strain}/{readorder}/{strain}.vcf", readorder = READORDER, allow_missing = True)
	output:
		expand("4_variant_call_comparisons/{depth}/{aligner}/pbsv/{strain}/{results_file}", results_file = COMPARISON_FILES, allow_missing = True)
	params:
		minsize="0",
		outdir="4_variant_call_comparisons/{depth}/{aligner}/pbsv/{strain}/"
	conda:  "yaml/compare_vcfs.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:30:00"
	shell:
		"python3 scripts/2_compare_vcfs/compare_vcfs.py {input} {params.outdir} pbsv"

# Compare the Sniffles SVs predicted from the original FASTQ files with the shuffled versions
rule compare_sniffles:
	input:
		expand("3_variant_calls/{depth}/{aligner}/sniffles/{strain}/{readorder}/{strain}.vcf", readorder = READORDER, allow_missing = True)
	output:
		expand("4_variant_call_comparisons/{depth}/{aligner}/sniffles/{strain}/{results_file}", results_file = COMPARISON_FILES, allow_missing = True)
	params:
		minsize="0",
		outdir="4_variant_call_comparisons/{depth}/{aligner}/sniffles/{strain}/"
	conda:  "yaml/compare_vcfs.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:30:00"
	shell:
		"python3 scripts/2_compare_vcfs/compare_vcfs.py {input} {params.outdir} sniffles"

# Compare the SVIM SVs predicted from the original FASTQ files with the shuffled versions
rule compare_svim:
	input:
		expand("3_variant_calls/{depth}/{aligner}/svim/{strain}/{readorder}/variants.vcf", readorder = READORDER, allow_missing = True)
	output:
		expand("4_variant_call_comparisons/{depth}/{aligner}/svim/{strain}/{results_file}", results_file = COMPARISON_FILES, allow_missing = True)
	params:
		quality = "0",
		minsize="0",
		outdir="4_variant_call_comparisons/{depth}/{aligner}/svim/{strain}/"
	conda:  "yaml/compare_vcfs.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:30:00"
	shell:
		"python3 scripts/2_compare_vcfs/compare_vcfs.py {input} {params.outdir} svim"

# Compare the pbsv SVs predicted from the original FASTQ files with the shuffled versions
rule compare_pbsv_original:
	input:
		original_vcf="3_variant_calls/{depth}/{aligner}/pbsv/{strain}/original/{strain}.vcf",
		shuffled_vcf="3_variant_calls/{depth}/{aligner}/pbsv/{strain}/{readorder}/{strain}.vcf"
	output:
		expand("4_variant_call_comparisons/{depth}/{aligner}/pbsv/{strain}/original_shuffled_comparisons/{readorder}/{results_file}", results_file = ORIGINAL_COMPARISON_FILES, allow_missing = True)
	params:
		quality = "0",
		minsize="0",
		outdir="4_variant_call_comparisons/{depth}/{aligner}/pbsv/{strain}/original_shuffled_comparisons/{readorder}/"
	conda:  "yaml/compare_vcfs.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/2_compare_vcfs/compare_original_shuffled.py {input.original_vcf} {input.shuffled_vcf} {params.outdir} pbsv"

# Compare the Sniffles SVs predicted from the original FASTQ files with the shuffled versions
rule compare_sniffles_original:
	input:
		original_vcf="3_variant_calls/{depth}/{aligner}/sniffles/{strain}/original/{strain}.vcf",
		shuffled_vcf="3_variant_calls/{depth}/{aligner}/sniffles/{strain}/{readorder}/{strain}.vcf"
	output:
		expand("4_variant_call_comparisons/{depth}/{aligner}/sniffles/{strain}/original_shuffled_comparisons/{readorder}/{results_file}", results_file = ORIGINAL_COMPARISON_FILES, allow_missing = True)
	params:
		quality = "0",
		minsize="0",
		outdir="4_variant_call_comparisons/{depth}/{aligner}/sniffles/{strain}/original_shuffled_comparisons/{readorder}/"
	conda:  "yaml/compare_vcfs.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/2_compare_vcfs/compare_original_shuffled.py {input.original_vcf} {input.shuffled_vcf} {params.outdir} sniffles"

# Compare the SVIM SVs predicted from the original FASTQ files with the shuffled versions
rule compare_svim_original:
	input:
		original_vcf="3_variant_calls/{depth}/{aligner}/svim/{strain}/original/variants.vcf",
		shuffled_vcf="3_variant_calls/{depth}/{aligner}/svim/{strain}/{readorder}/variants.vcf"
	output:
		expand("4_variant_call_comparisons/{depth}/{aligner}/svim/{strain}/original_shuffled_comparisons/{readorder}/{results_file}", results_file = ORIGINAL_COMPARISON_FILES, allow_missing = True)
	params:
		quality = "0",
		minsize="0",
		outdir="4_variant_call_comparisons/{depth}/{aligner}/svim/{strain}/original_shuffled_comparisons/{readorder}/"
	conda:  "yaml/compare_vcfs.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"python3 scripts/2_compare_vcfs/compare_original_shuffled.py {input.original_vcf} {input.shuffled_vcf} {params.outdir} svim"

# Summarize the SV call agreement per strain
rule summarize_vcf_agreement_per_strain:
	input:
		pbsv_files=expand("4_variant_call_comparisons/{depth}/{aligner}/pbsv/{strain}/original_shuffled_comparisons/{shuffle_dir}/results_proportions.csv", shuffle_dir = SHUFFLED_DIRS, allow_missing = True),
		sniffles_files=expand("4_variant_call_comparisons/{depth}/{aligner}/sniffles/{strain}/original_shuffled_comparisons/{shuffle_dir}/results_proportions.csv", shuffle_dir = SHUFFLED_DIRS, allow_missing = True),
		svim_files=expand("4_variant_call_comparisons/{depth}/{aligner}/svim/{strain}/original_shuffled_comparisons/{shuffle_dir}/results_proportions.csv", shuffle_dir = SHUFFLED_DIRS, allow_missing = True),
	output:
		"5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/{strain}_results.csv"
	params:
		outdir="5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/"
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:20:00"
	shell:
		"python scripts/3_summarize_vcf_comparisons/summarize_vcf_agreement_per_strain.py -p {input.pbsv_files} -s {input.sniffles_files} -v {input.svim_files} -o {params.outdir} -n {wildcards.strain}"

# Summarize the comparisons between the original and shuffled vcf files for each strain without including metadata
rule summarize_vcf_agreement_per_strain_no_meta:
	input:
		pbsv_files=expand("4_variant_call_comparisons/{depth}/{aligner}/pbsv/{strain}/original_shuffled_comparisons/{shuffle_dir}/results_proportions_no_meta.csv", shuffle_dir = SHUFFLED_DIRS, allow_missing = True),
		sniffles_files=expand("4_variant_call_comparisons/{depth}/{aligner}/sniffles/{strain}/original_shuffled_comparisons/{shuffle_dir}/results_proportions_no_meta.csv", shuffle_dir = SHUFFLED_DIRS, allow_missing = True),
		svim_files=expand("4_variant_call_comparisons/{depth}/{aligner}/svim/{strain}/original_shuffled_comparisons/{shuffle_dir}/results_proportions_no_meta.csv", shuffle_dir = SHUFFLED_DIRS, allow_missing = True),
	output:
		"5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/{strain}_results_no_meta.csv"
	params:
		outdir="5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/"
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:20:00"
	shell:
		"python scripts/3_summarize_vcf_comparisons/summarize_vcf_agreement_per_strain.py -p {input.pbsv_files} -s {input.sniffles_files} -v {input.svim_files} -o {params.outdir} -n {wildcards.strain} --no_meta"

# Contactenate vcf agreement stats per strain
rule concatenate_per_strain_vcf_agreement:
	input:
		expand("5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/{strain}_results.csv", strain = ALL_STRAINS, allow_missing = True)
	output:
		"5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/results_combined.csv"
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:20:00"
	shell:
		"python scripts/3_summarize_vcf_comparisons/concatenate_agreement_per_strain_results.py -i {input} -o {output}"

# Contactenate vcf agreement stats per strain for xyz.no_meta.csv files
rule concatenate_per_strain_vcf_agreement_no_meta:
	input:
		expand("5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/{strain}_results_no_meta.csv", strain = ALL_STRAINS, allow_missing = True)
	output:
		"5_vcf_comparison_summaries/{depth}/{aligner}/original_shuffled_agreement_per_strain/results_combined_no_meta.csv"
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:20:00"
	shell:
		"python scripts/3_summarize_vcf_comparisons/concatenate_agreement_per_strain_results.py -i {input} -o {output}"

# Concatenate agreement stats for all strains and replicates
rule concatenate_vcf_agreement_totals:
	input:
		expand("4_variant_call_comparisons/{depth}/{aligner}/{caller}/{strain}/original_shuffled_comparisons/{readorder}/results_proportions.csv", strain = ALL_STRAINS, readorder = SHUFFLED_DIRS, allow_missing = True)
	output:
		"5_vcf_comparison_summaries/{depth}/{aligner}/{caller}/csv_totals_all_strains_replicates_jaccard.csv",
		"5_vcf_comparison_summaries/{depth}/{aligner}/{caller}/csv_totals_all_strains_replicates_differences.csv",
	params:
		input_format = "4_variant_call_comparisons/{depth}/{aligner}/{caller}/STRAIN/original_shuffled_comparisons/shuffledREPLICATE/results_proportions.csv",
		outdir="5_vcf_comparison_summaries/{depth}/{aligner}/{caller}/",
		replicates = "5"
	conda:  "yaml/pandas.yaml"
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="00:30:00"
	shell:
		"python3 scripts/3_summarize_vcf_comparisons/summarize_original_shuffled_comparisons.py {params.input_format} {params.outdir}"
