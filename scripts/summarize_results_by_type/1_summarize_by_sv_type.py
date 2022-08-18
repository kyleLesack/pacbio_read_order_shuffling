#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:21:34 2022

@author: kyle
"""
import os
from collections import defaultdict
from itertools import chain
import statistics

#os.chdir('/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling')

SNIFFLES_PATH = "3_variant_calls/ALIGNER/sniffles/STRAIN/REP/summary/"
#SVIM_PATH_QUAL_0 = "3_variant_calls/ALIGNER/svim/STRAIN/REP/QUAL_0/summary/"
SVIM_PATH_QUAL_15 = "3_variant_calls/ALIGNER/svim/STRAIN/REP/QUAL_15/summary/"
PBSV_PATH = "3_variant_calls/ALIGNER/pbsv/STRAIN/REP/summary/"

ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]
SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["pbmm2"]
PBSV_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV"]
SNIFFLES_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV"]
SVIM_SV_TYPES = ["BND", "DEL", "DUP:INT","DUP:TANDEM", "INS","INV"]

SUBSAMPLE_DEPTHS = ["10X","20X","40X","60X"]

SV_Intersection_HEADER = "SV TYPE,Intersection,Original only,Shuffled only"
BREAKPOINT_HEADER = "SV TYPE,Intersection,Same breakpoints,Different breakpoints"

def import_lines(filename, aligner,strain, caller_sv_types):
	with open(filename) as f:
		my_lines = f.readlines()
		
	sv_type_dict = defaultdict(list) # Dictionary to store the agreement stats for each sv type. 
	for line in my_lines[1:]: # Skip header
		split_lines = line.strip().split(",")
		sv_type = split_lines[1]
		new_line = [split_lines[0]]
		for x in split_lines[2:]:
			new_line.append(x)
		if sv_type in caller_sv_types:
			sv_type_dict[sv_type].append(new_line)

	# Dictionary to store the agreement stats for each sv type. If a given SV type has multiple filters, the agreements values are summed together.
	sv_type_dict_totals = defaultdict(list) 
	for x in sv_type_dict.keys():
		sv_agreement_lines = sv_type_dict[x]
		line_count = len(sv_agreement_lines) # Number of values for the key. There are > 1 if multiple filter types
		if line_count > 1:
			Intersection_total = 0
			same_breakpoint_total = 0
			different_breakpoint_total = 0
			total_filters = []
			for line in sv_agreement_lines:
				total_filters.append(line[0])
				Intersection_total += int(line[1])
				same_breakpoint_total += int(line[2])
				different_breakpoint_total += int(line[3])
			total_filters_joined = ";".join(total_filters)
			new_line = [strain, str(total_filters_joined),Intersection_total, same_breakpoint_total,different_breakpoint_total]
			sv_type_dict_totals[x].append(new_line)
			#print(new_line)
		else:
			sv_agreement_line = list(chain.from_iterable(sv_agreement_lines))
			sv_agreement_line.insert(0,strain)
			sv_type_dict_totals[x].append(sv_agreement_line)
	
	return(sv_type_dict_totals)

# Function calculates the agreement statistics for each caller/aligner combo. 
# Parameters: caller = caller name ,caller_csv_path = path to caller csv files, aligners = aligners used to generate alignments for the caller, caller_sv_types = sv types predicted by the caller, replicate = the shuffled genome to compare to the original
def get_agreement_summary(caller,caller_csv_path, aligners, caller_sv_types, replicate, inputfile):
	# caller_aligner_dict_bp: dictionary to store the breakpoint agreement data imported from csv files
	# caller_aligner_dict_bp: keys = caller/aligner combo and the svtype
	# caller_aligner_dict_bp: values = filters, Intersection, same breakpoints, different breakpoints
	caller_aligner_dict_bp = defaultdict(lambda: defaultdict(list)) 
	
	# Loop is used to populate caller_aligner_dict_bp
	for aligner in aligners: # analyze results for each different aligner used for the caller
		caller_aligner_combo = caller + "-" + aligner
		for strain in ALL_STRAINS:
			caller_csv_file = caller_csv_path.replace("ALIGNER",aligner).replace("STRAIN",strain).replace("REP",str(replicate)) + inputfile
			if os.path.isfile(caller_csv_file):
				# breakpoint_agreement_dict stores csv file information
				# breakpoint_agreement_dict: keys = Variant type, values = filters, Intersection, same breakpoints, different breakpoints
				breakpoint_agreement_dict = import_lines(caller_csv_file, aligner,strain, caller_sv_types) 

				# The following loop is used to populate the caller_aligner_dict_bp dictionary with the breakpoint_agreement_dict data
				for svtype in breakpoint_agreement_dict.keys():
					#print(svtype)
					# get the filters, Intersection, same breakpoints, different breakpoints for a given sample and store in list
					caller_aligner_bp_line = list(chain.from_iterable(breakpoint_agreement_dict[svtype])) 
					caller_aligner_dict_bp[caller_aligner_combo][svtype].append(caller_aligner_bp_line)
			else:
				print("File not found: " + caller_csv_file)
	
	# caller_aligner_dict_bp_means stores the mean agreement statistic values for each "caller/aligner and svtype" combination
	# caller_aligner_dict_bp_means: keys = caller/aligner combination and svtype
	# caller_aligner_dict_bp_means: values = mean counts for the Intersection, calls with breakpoint agreement, calls with breakpoint disagreement
	caller_aligner_dict_bp_means = defaultdict(lambda: defaultdict(list)) # Dictionary to store the mean 
	#caller_aligner_dict_bp_means = defaultdict(lambda: defaultdict)
	
	for caller_aligner in caller_aligner_dict_bp.keys(): # Loop to calculate means for each statistic
		for svtype in caller_aligner_dict_bp[caller_aligner]:
			Intersection_values = []
			same_breakpoint_values = []
			different_breakpoint_values = []
			filters = set()

			for stat_line in caller_aligner_dict_bp[caller_aligner][svtype]:
				filters.add(stat_line[1])
				Intersection_values.append(int(stat_line[2]))
				same_breakpoint_values.append(int(stat_line[3]))
				different_breakpoint_values.append(int(stat_line[4]))

			mean_intersection = statistics.mean(Intersection_values)
			mean_same_breakpoints = statistics.mean(same_breakpoint_values)
			mean_different_breakpoints = statistics.mean(different_breakpoint_values)
			summary_line =f'{mean_intersection},{mean_same_breakpoints},{mean_different_breakpoints}'
			#summary_line =f'{mean_intersection:.2f},{mean_same_breakpoints:.2f},{mean_different_breakpoints:.2f}'
			#summary_line =f'{mean_intersection:.0f},{mean_same_breakpoints:.0f},{mean_different_breakpoints:.0f}'
			caller_aligner_dict_bp_means[caller_aligner][svtype]= summary_line

	return(caller_aligner_dict_bp_means)


# Write results to disk
def write_results(result_dict, outputpath, header):
	if not os.path.exists(outputpath):
		os.makedirs(outputpath)
	for caller_aligner in result_dict.keys():
		outfile = outputpath + caller_aligner + ".csv"
		with open(outfile, 'w', newline='\n') as f:
			f.write(header+"\n")
			for svtype in result_dict[caller_aligner]:
				csv_line = svtype + "," + result_dict[caller_aligner][svtype] + "\n"
				f.write(csv_line)

# Summarize breakpoint agreement for all strains
sniffles_breakpoint_means = get_agreement_summary("sniffles",SNIFFLES_PATH, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_breakpoints.csv")
#svim_qual0_breakpoint_means = get_agreement_summary("svim", SVIM_PATH_QUAL_0, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_breakpoints.csv")	
svim_qual15_breakpoint_means = get_agreement_summary("svim", SVIM_PATH_QUAL_15, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_breakpoints.csv")
pbsv_breakpoint_means = get_agreement_summary("pbsv",PBSV_PATH, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_breakpoints.csv")

write_results(sniffles_breakpoint_means, "4_results/full_depth/breakpoint_agreement/sniffles/", BREAKPOINT_HEADER)
#write_results(svim_qual0_breakpoint_means, "4_results/full_depth/breakpoint_agreement/svim/qual_0/", BREAKPOINT_HEADER)
write_results(svim_qual15_breakpoint_means, "4_results/full_depth/breakpoint_agreement/svim/qual_15/", BREAKPOINT_HEADER)
write_results(pbsv_breakpoint_means, "4_results/full_depth/breakpoint_agreement/pbsv/", BREAKPOINT_HEADER)

# Subsampled Analysis - Summarize breakpoint agreement for all strains
for depth in SUBSAMPLE_DEPTHS:
	sniffles_subsampled_path = SNIFFLES_PATH
	subdepth = "/subsampled/" + depth + "/"
	sniffles_subsampled_path = sniffles_subsampled_path.replace("/", subdepth, 1)
	sniffles_breakpoint_means = get_agreement_summary("sniffles",sniffles_subsampled_path, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_breakpoints.csv")
	outdir = "4_results" + subdepth + "breakpoint_agreement/sniffles/"
	write_results(sniffles_breakpoint_means, outdir, BREAKPOINT_HEADER)
	
	pbsv_subsampled_path = PBSV_PATH
	subdepth = "/subsampled/" + depth + "/"
	pbsv_subsampled_path = pbsv_subsampled_path.replace("/", subdepth, 1)
	pbsv_breakpoint_means = get_agreement_summary("pbsv",pbsv_subsampled_path, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_breakpoints.csv")
	outdir = "4_results" + subdepth + "breakpoint_agreement/pbsv/"
	write_results(pbsv_breakpoint_means, outdir, BREAKPOINT_HEADER)
	
	#svim_subsampled_path = SVIM_PATH_QUAL_0
	#subdepth = "/subsampled/" + depth + "/"
	#svim_subsampled_path = svim_subsampled_path.replace("/", subdepth, 1)
	#svim_breakpoint_means = get_agreement_summary("svim",svim_subsampled_path, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_breakpoints.csv")
	#outdir = "4_results" + subdepth + "breakpoint_agreement/svim/qual_0/"
	#write_results(svim_breakpoint_means, outdir, BREAKPOINT_HEADER)
	
	svim_subsampled_path = SVIM_PATH_QUAL_15
	subdepth = "/subsampled/" + depth + "/"
	svim_subsampled_path = svim_subsampled_path.replace("/", subdepth, 1)
	svim_breakpoint_means = get_agreement_summary("svim",svim_subsampled_path, SNIFFLES_SVIM_ALIGNERS , SVIM_SV_TYPES, 1, "summary_breakpoints.csv")
	outdir = "4_results" + subdepth + "breakpoint_agreement/svim/qual_15/"
	write_results(svim_breakpoint_means, outdir, BREAKPOINT_HEADER)	

# Summarize structural Intersection agreement for all strains
sniffles_sv_Intersection_means = get_agreement_summary("sniffles",SNIFFLES_PATH, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_total_svs.csv")
#svim_qual0_sv_Intersection_means = get_agreement_summary("svim", SVIM_PATH_QUAL_0, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_total_svs.csv")	
svim_qual15_sv_Intersection_means = get_agreement_summary("svim", SVIM_PATH_QUAL_15, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_total_svs.csv")
pbsv_sv_Intersection_means = get_agreement_summary("pbsv",PBSV_PATH, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_total_svs.csv")

write_results(sniffles_sv_Intersection_means, "4_results/full_depth/sv_intersection_agreement/sniffles/", SV_Intersection_HEADER)
#write_results(svim_qual0_sv_Intersection_means, "4_results/full_depth/sv_intersection_agreement/svim/qual_0/", SV_Intersection_HEADER)
write_results(svim_qual15_sv_Intersection_means, "4_results/full_depth/sv_intersection_agreement/svim/qual_15/", SV_Intersection_HEADER)
write_results(pbsv_sv_Intersection_means, "4_results/full_depth/sv_intersection_agreement/pbsv/", SV_Intersection_HEADER)

# Subsampled Analysis - Summarize structural variant Intersection agreement for all strains
for depth in SUBSAMPLE_DEPTHS:
	sniffles_subsampled_path = SNIFFLES_PATH
	subdepth = "/subsampled/" + depth + "/"
	sniffles_subsampled_path = sniffles_subsampled_path.replace("/", subdepth, 1)
	sniffles_sv_Intersection_means = get_agreement_summary("sniffles",sniffles_subsampled_path, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_total_svs.csv")
	outdir = "4_results" + subdepth + "sv_intersection_agreement/sniffles/"
	write_results(sniffles_sv_Intersection_means, outdir, SV_Intersection_HEADER)
	
	pbsv_subsampled_path = PBSV_PATH
	subdepth = "/subsampled/" + depth + "/"
	pbsv_subsampled_path = pbsv_subsampled_path.replace("/", subdepth, 1)
	pbsv_sv_Intersection_means = get_agreement_summary("pbsv",pbsv_subsampled_path, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_total_svs.csv")
	outdir = "4_results" + subdepth + "sv_intersection_agreement/pbsv/"
	write_results(pbsv_sv_Intersection_means, outdir, SV_Intersection_HEADER)
	
	#svim_subsampled_path = SVIM_PATH_QUAL_0
	#subdepth = "/subsampled/" + depth + "/"
	#svim_subsampled_path = svim_subsampled_path.replace("/", subdepth, 1)
	#svim_sv_Intersection_means = get_agreement_summary("svim",svim_subsampled_path, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_total_svs.csv")
	#outdir = "4_results" + subdepth + "sv_intersection_agreement/svim/qual_0/"
	#write_results(svim_sv_Intersection_means, outdir, SV_Intersection_HEADER)
	
	svim_subsampled_path = SVIM_PATH_QUAL_15
	subdepth = "/subsampled/" + depth + "/"
	svim_subsampled_path = svim_subsampled_path.replace("/", subdepth, 1)
	svim_sv_Intersection_means = get_agreement_summary("svim",svim_subsampled_path, SNIFFLES_SVIM_ALIGNERS , SVIM_SV_TYPES, 1, "summary_total_svs.csv")
	outdir = "4_results" + subdepth + "sv_intersection_agreement/svim/qual_15/"
	write_results(svim_sv_Intersection_means, outdir, SV_Intersection_HEADER)	
