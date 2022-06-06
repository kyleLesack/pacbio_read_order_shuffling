include: "0_input/includes/Snakefile.no_shuffle"

BIONET_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147"]
SPLIT_STRAINS = ["DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"] # Large strains that took too long to shuffle. I split the Fastq files before shuffling and concatenation.
NO_SPLIT_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768"]
ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]
REFERENCE = "/work/wasmuth_lab/mrkyle/sv_calling_pipeline2/1_prepare_reference/output/c_elegans.PRJNA13758.WS263/c_elegans.PRJNA13758.WS263.genomic.fa"
REFERENCE_SAW = "/work/wasmuth_lab/mrkyle/sv_calling_pipeline2/1_prepare_reference/output/c_elegans.PRJNA13758.WS263.sawriter/c_elegans.PRJNA13758.WS263.genomic.fa"
SAM_ALIGNERS = ["ngmlr","minimap2"] # Aligners that output sam files
SPLIT_DIRS = ["no_split", "needs_split"]
SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
SNIFFLES_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["pbmm2"]
ALL_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]


NGMLRDICT10X = {"DL238": "0.06", "DRR142768": "0.13", "ECA36": "0.05", "ECA396": "0.07", "EG4725": "0.05", "JU1400": "0.09", "JU2526": "0.08", "JU2600": "0.06", "JU310": "0.06", "MY2147": "0.05", "MY2693": "0.08", "NIC2": "0.08", "NIC526": "0.07", "QX1794": "0.08","XZ1516": "0.13"}
MINIMAP2XDICT10X = {'DL238': '0.05', 'DRR142768': '0.12', 'ECA36': '0.04', 'ECA396': '0.07', 'EG4725': '0.04', 'JU1400': '0.08', 'JU2526': '0.08', 'JU2600': '0.05', 'JU310': '0.05', 'MY2147': '0.05', 'MY2693': '0.07', 'NIC2': '0.08', 'NIC526': '0.06', 'QX1794': '0.07','XZ1516': '0.12'}
PBMM2XDICT10X = {'DL238': '0.06', 'DRR142768': '0.12', 'ECA36': '0.05', 'ECA396': '0.07', 'EG4725': '0.04', 'JU1400': '0.08', 'JU2526': '0.08', 'JU2600': '0.05', 'JU310': '0.06', 'MY2147': '0.05', 'MY2693': '0.07', 'NIC2': '0.08', 'NIC526': '0.06', 'QX1794': '0.08','XZ1516': '0.13'}


rule all:
	input:
		expand("1_fq_processing/all_strains/shuffled/{no_split_strain}/{replicate}/{no_split_strain}_shuffled.fastq", no_split_strain=NO_SPLIT_STRAINS, replicate=[1,2,3,4,5]),
		expand("1_fq_processing/all_strains/shuffled/{split_strain}/{replicate}/{split_strain}_shuffled.fastq", split_strain=SPLIT_STRAINS, replicate=[1,2,3,4,5]),
		expand("2_alignments/{aligner}/{strain}/ORIGINAL/{strain}.sam", aligner=SAM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam", aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/{aligner}/{strain}/ORIGINAL/{strain}_sorted.bam.bai", aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/{aligner}/{strain}/ORIGINAL/{strain}_depth.txt", aligner=SAM_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/subsampled/10X/ngmlr/{strain}/ORIGINAL/{strain}_sorted.bam",aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_shuffled.sam", aligner=SAM_ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam", aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_shuffled_sorted.bam.bai", aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("2_alignments/{aligner}/{strain}/{replicate}/{strain}_depth.txt", aligner=SAM_ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("2_alignments/subsampled/10X/ngmlr/{strain}/{replicate}/{strain}_shuffled_sorted.bam",aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS, replicate=[1]),
		#expand("2_alignments/pbmm2/{strain}/{replicate}/{strain}_shuffled.bam", strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		#expand("2_alignments/pbmm2/{strain}/{replicate}/{strain}_depth.txt", strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf", aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/svim/{strain}/{replicate}/variants.vcf", aligner=SVIM_ALIGNERS, strain=ALL_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf", aligner=PBSV_ALIGNERS, strain=ALL_STRAINS),
		expand("3_variant_calls/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf", aligner=SNIFFLES_ALIGNERS, strain=ALL_STRAINS),
		expand("3_variant_calls/{aligner}/svim/{strain}/ORIGINAL/variants.vcf", aligner=SVIM_ALIGNERS, strain=ALL_STRAINS),
		expand("3_variant_calls/{aligner}/svim/{strain}/{replicate}/QUAL_0/summary/summary_total_svs.csv", aligner=SVIM_ALIGNERS, strain=BIONET_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/svim/{strain}/{replicate}/QUAL_0/summary/summary_breakpoints.csv", aligner=SVIM_ALIGNERS, strain=BIONET_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv", aligner=SNIFFLES_ALIGNERS, strain=BIONET_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv", aligner=SNIFFLES_ALIGNERS, strain=BIONET_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv", aligner=PBSV_ALIGNERS, strain=BIONET_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv", aligner=PBSV_ALIGNERS, strain=BIONET_STRAINS, replicate=[1,2,3,4,5]),
		expand("3_variant_calls/{aligner}/pbsv/DL238/1/DL238.vcf", aligner=ALL_ALIGNERS)		
		


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

rule pbmm2:
	input:
		"1_fq_processing/all_strains/shuffled/{strain}/{REP}/{strain}_shuffled.fastq"
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


rule subsamplengmlrto10x:
	input:
		"2_alignments/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	output:
		"2_alignments/subsampled/10X/ngmlr/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	params:
		value=lambda wcs: NGMLRDICT10X[wcs.strain]
	#        wildcard_constraints:
	#               aligner = "|".join(SAM_ALIGNERS)
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} -b {input} > {output}"""

rule pbsvdiscover:
	input:
		"2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled.bam"
	output:
		temp("3_variant_calls/{aligner}/{strain}/{REP}/{strain}.svsig.gz")
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
		"3_variant_calls/{aligner}/{strain}/{REP}/{strain}.svsig.gz"
	output:
		"3_variant_calls/{aligner}/pbsv/{strain}/{REP}/{strain}.vcf"
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
		"3_variant_calls/{aligner}/svim/{strain}/{REP}/variants.vcf"
	wildcard_constraints:
		aligner = "|".join(SVIM_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	params:
		outdir="3_variant_calls/{aligner}/svim/{strain}/{REP}/"
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
		aligner = "|".join(SNIFFLES_ALIGNERS)
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
		"3_variant_calls/{aligner}/sniffles/{strain}/{REP}/{strain}.vcf"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_ALIGNERS),
		REP = "|".join(["1","2","3","4","5"])
	params:
		bamfile="2_alignments/{aligner}/{strain}/{REP}/{strain}_shuffled_sorted.bam"
	conda:  "yaml/sniffles2.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="12:00:00"
	shell:
		"sniffles --input {params.bamfile} -t 8 --vcf {output}"

rule compare_svim:
	input:
		"3_variant_calls/{aligner}/svim/{strain}/{replicate}/variants.vcf",
		"3_variant_calls/{aligner}/svim/{strain}/ORIGINAL/variants.vcf"
	output:
		"3_variant_calls/{aligner}/svim/{strain}/{replicate}/QUAL_0/summary/summary_total_svs.csv",
		"3_variant_calls/{aligner}/svim/{strain}/{replicate}/QUAL_0/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(SVIM_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
	params:
		originalnum = "3_variant_calls/{aligner}/svim/{strain}/ORIGINAL/variants.vcf",
		shuffled = "3_variant_calls/{aligner}/svim/{strain}/{replicate}/variants.vcf"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
			"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.originalnum} svim --minsize 0 --min_qual_svim 0"

rule compare_sniffles:
	input:
		"3_variant_calls/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf",
		"3_variant_calls/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf"
	output:
		"3_variant_calls/{aligner}/sniffles/{strain}/{replicate}/summary/summary_total_svs.csv",
		"3_variant_calls/{aligner}/sniffles/{strain}/{replicate}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(SNIFFLES_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
	params:
		originalnum = "3_variant_calls/{aligner}/sniffles/{strain}/ORIGINAL/{strain}.vcf",
		shuffled = "3_variant_calls/{aligner}/sniffles/{strain}/{replicate}/{strain}.vcf"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
			"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.originalnum} sniffles --minsize 0"

rule compare_pbsv:
	input:
		"3_variant_calls/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf",
		"3_variant_calls/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf"
	output:
		"3_variant_calls/{aligner}/pbsv/{strain}/{replicate}/summary/summary_total_svs.csv",
		"3_variant_calls/{aligner}/pbsv/{strain}/{replicate}/summary/summary_breakpoints.csv"
	wildcard_constraints:
		aligner = "|".join(PBSV_ALIGNERS),
		replicate = "|".join(["1","2","3","4","5"]),
	params:
		originalnum = "3_variant_calls/{aligner}/pbsv/{strain}/ORIGINAL/{strain}.vcf",
		shuffled = "3_variant_calls/{aligner}/pbsv/{strain}/{replicate}/{strain}.vcf"
	conda:  "yaml/pybedtools.yaml"
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="01:00:00"
	shell:
			"python3 scripts/compare_shuffled_2_original.py {params.shuffled} {params.originalnum} pbsv --minsize 0"

rule pbsvdiscover_ngmlr:
	input:
		"2_alignments/ngmlr/DL238/1/DL238_shuffled_sorted.bam"
	output:
		temp("3_variant_calls/ngmlr/pbsv/DL238/1/DL238.svsig.gz")
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="02:00:00"
	shell:
		"pbsv discover -s {input} {output}"

rule pbsvcall_ngmlr:
	input:
		"3_variant_calls/ngmlr/pbsv/DL238/1/DL238.svsig.gz"
	output:
		"3_variant_calls/ngmlr/pbsv/DL238/1/DL238.vcf"
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),   
		time="08:00:00"
	shell:
		"pbsv call -j 8 {REFERENCE_SAW} {input} {output}"

rule pbsvdiscover_minimap2:
	input:
		"2_alignments/minimap2/DL238/1/DL238_shuffled_sorted.bam"
	output:
		temp("3_variant_calls/minimap2/pbsv/DL238/1/DL238.svsig.gz")
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),
		time="02:00:00"
	shell:
		"pbsv discover -s {input} {output}"

rule pbsvcall_minimap2:
	input:
		"3_variant_calls/minimap2/pbsv/DL238/1/DL238.svsig.gz"
	output:
		"3_variant_calls/minimap2/pbsv/DL238/1/DL238.vcf"
	conda:  "yaml/pbsv.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 50000 + ((attempt - 1) * 10000),   
		time="08:00:00"
	shell:
		"pbsv call -j 8 {REFERENCE_SAW} {input} {output}"

