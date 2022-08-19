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

os.chdir('/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling')

SNIFFLES_PATH = "3_variant_calls/DEPTH/ALIGNER/sniffles/STRAIN/REP/summary/"
SVIM_PATH_QUAL_15 = "3_variant_calls/DEPTH/ALIGNER/svim/STRAIN/REP/QUAL_15/summary/"
PBSV_PATH = "3_variant_calls/DEPTH/ALIGNER/pbsv/STRAIN/REP/summary/"

ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DRR142768", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]
SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["pbmm2"]
PBSV_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV","ALL_SVS"]
SNIFFLES_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV","ALL_SVS"]
SVIM_SV_TYPES = ["BND", "DEL", "DUP:INT","DUP:TANDEM", "INS","INV","ALL_SVS"]

SUBSAMPLE_DEPTHS = ["10X","20X","40X","60X"]

SV_INTERSECTION_HEADER = "SV TYPE,Intersection,Original only,Shuffled only,Unique,Unique proportion"
BREAKPOINT_HEADER = "SV TYPE,Intersection,Same breakpoints,Different breakpoints,Discordant proportion"

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
		else:
			print("SV Type excluded: " + sv_type)
	# Dictionary to store the agreement stats for each sv type. If a given SV type has multiple filters, the agreements values are summed together.
	sv_type_dict_totals = defaultdict(list) 
	for x in sv_type_dict.keys():
		sv_call_lines = sv_type_dict[x]
		line_count = len(sv_call_lines) # Number of values for the key. There are > 1 if multiple filter types
		
		if "summary_breakpoints" in filename:
			if line_count > 1:
				intersection_total = 0
				same_breakpoint_total = 0
				different_breakpoint_total = 0
				total_filters = []
				for line in sv_call_lines:
					total_filters.append(line[0])
					intersection_total += int(line[1])
					same_breakpoint_total += int(line[2])
					different_breakpoint_total += int(line[3])
				discordant_proportion = different_breakpoint_total/intersection_total 
				total_filters_joined = ";".join(total_filters)
				new_line = [strain, str(total_filters_joined),intersection_total, same_breakpoint_total,different_breakpoint_total,discordant_proportion]
				sv_type_dict_totals[x].append(new_line)
				
			else:
				sv_agreement_line = list(chain.from_iterable(sv_call_lines))
				sv_agreement_line.insert(0,strain)
				sv_type_dict_totals[x].append(sv_agreement_line)

		elif "summary_total_svs" in filename:
			if line_count > 1:
				intersection_total = 0
				original_only_total = 0
				shuffled_only_total = 0
				unique_total = 0
				total_filters = []
				for line in sv_call_lines:
					total_filters.append(line[0])
					intersection_total += int(line[1])
					original_only_total += int(line[2])
					shuffled_only_total += int(line[3])
					unique_total += int(line[4])
				total_calls = intersection_total + unique_total
				unique_proportion = unique_total / total_calls
				total_filters_joined = ";".join(total_filters)
				new_line = [strain, str(total_filters_joined),intersection_total, original_only_total,shuffled_only_total,unique_total, unique_proportion]
				sv_type_dict_totals[x].append(new_line)
	
			else:
				sv_agreement_line = list(chain.from_iterable(sv_call_lines))
				sv_agreement_line.insert(0,strain)
				sv_type_dict_totals[x].append(sv_agreement_line)
	return(sv_type_dict_totals)

# Function calculates the agreement statistics for each caller/aligner combo. 
# Parameters: caller = caller name ,caller_csv_path = path to caller csv files, aligners = aligners used to generate alignments for the caller, caller_sv_types = sv types predicted by the caller, replicate = the shuffled genome to compare to the original
def get_agreement_summary(caller,caller_csv_path, aligners, caller_sv_types, replicate, inputfile):
	#print(inputfile)
	# caller_aligner_dict_bp: dictionary to store the sv overlap or breakpoint agreement data imported from csv files
	# caller_aligner_dict_bp: keys = caller/aligner combo and the svtype
	# caller_aligner_dict_bp: values = filters, intersection, same breakpoints, different breakpoints
	caller_aligner_dict_bp = defaultdict(lambda: defaultdict(list)) 
	
	# Loop is used to populate caller_aligner_dict_bp
	for aligner in aligners: # analyze results for each different aligner used for the caller
		caller_aligner_combo = caller + "-" + aligner
		for strain in ALL_STRAINS:
			caller_csv_file = caller_csv_path.replace("ALIGNER",aligner).replace("STRAIN",strain).replace("REP",str(replicate)) + inputfile

			if os.path.isfile(caller_csv_file):
				# sv_agreement_dict stores csv file information
				# sv_agreement_dict: keys = Variant type, values = filters, intersection, same breakpoints or original only, different breakpoints or shuffled only, discordant proportion or total unique, unique proportion (only if summary_total_svs.csv was the csv file)
				sv_agreement_dict = import_lines(caller_csv_file, aligner,strain, caller_sv_types) 
				# The following loop is used to populate the caller_aligner_dict_bp dictionary with the sv_agreement_dict data
				for svtype in sv_agreement_dict.keys():
					caller_aligner_bp_line = list(chain.from_iterable(sv_agreement_dict[svtype])) 
					caller_aligner_dict_bp[caller_aligner_combo][svtype].append(caller_aligner_bp_line)
			else:
				print("File not found: " + caller_csv_file)
	
	# caller_aligner_dict_bp_means stores the mean agreement statistic values for each "caller/aligner and svtype" combination
	# caller_aligner_dict_bp_means: keys = caller/aligner combination and svtype
	# caller_aligner_dict_bp_means: values = mean counts for the sv_agreement_dict fields
	caller_aligner_dict_bp_means = defaultdict(lambda: defaultdict(list)) # Dictionary to store the mean 
	#caller_aligner_dict_bp_means = defaultdict(lambda: defaultdict)
	
	for caller_aligner in caller_aligner_dict_bp.keys(): # Loop to calculate means for each statistic
		if "summary_total_svs.csv" in inputfile:
			
			for svtype in caller_aligner_dict_bp[caller_aligner]:
				#print(svtype)
				intersection_values = []
				original_only_values = []
				shuffled_only_values = []
				unique_values = []
				unique_proportion_values = []
				filters = set()
	
				for stat_line in caller_aligner_dict_bp[caller_aligner][svtype]:
					#print(stat_line)
					filters.add(stat_line[1])
					intersection_values.append(int(stat_line[2]))
					original_only_values.append(int(stat_line[3]))
					shuffled_only_values.append(int(stat_line[4]))
					unique_values.append(float(stat_line[5]))
					unique_proportion_values.append(float(stat_line[6]))
	
				mean_intersection = statistics.mean(intersection_values)
				mean_original_only = statistics.mean(original_only_values)
				mean_shuffled_only = statistics.mean(shuffled_only_values)
				mean_unique = statistics.mean(unique_values)
				mean_unique_proportion = statistics.mean(unique_proportion_values)
				summary_line =f'{mean_intersection},{mean_original_only},{mean_shuffled_only},{mean_unique},{mean_unique_proportion}'
				caller_aligner_dict_bp_means[caller_aligner][svtype]= summary_line

		elif "summary_breakpoints.csv" in inputfile:
			for svtype in caller_aligner_dict_bp[caller_aligner]:
				intersection_values = []
				same_breakpoint_values = []
				different_breakpoint_values = []
				discordant_proportion_values = []
				filters = set()
	
				for stat_line in caller_aligner_dict_bp[caller_aligner][svtype]:
					filters.add(stat_line[1])
					intersection_values.append(int(stat_line[2]))
					same_breakpoint_values.append(int(stat_line[3]))
					different_breakpoint_values.append(int(stat_line[4]))
					discordant_proportion_values.append(float(stat_line[5]))
	
				mean_intersection = statistics.mean(intersection_values)
				mean_same_breakpoints = statistics.mean(same_breakpoint_values)
				mean_different_breakpoints = statistics.mean(different_breakpoint_values)
				mean_discordant_proportion = statistics.mean(discordant_proportion_values)
				summary_line =f'{mean_intersection},{mean_same_breakpoints},{mean_different_breakpoints},{mean_discordant_proportion}'
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

# Summarize structural variant agreement for all strains
sniffles_path = SNIFFLES_PATH.replace("DEPTH","full_depth")
sniffles_sv_intersection_means = get_agreement_summary("sniffles",sniffles_path, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_total_svs.csv")
svim_path = SVIM_PATH_QUAL_15.replace("DEPTH","full_depth")
svim_qual15_sv_intersection_means = get_agreement_summary("svim", svim_path, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_total_svs.csv")
pbsv_path = PBSV_PATH.replace("DEPTH","full_depth")
pbsv_sv_intersection_means = get_agreement_summary("pbsv",pbsv_path, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_total_svs.csv")

write_results(sniffles_sv_intersection_means, "4_results/full_depth/sv_intersection_agreement/sniffles/", SV_INTERSECTION_HEADER)
write_results(svim_qual15_sv_intersection_means, "4_results/full_depth/sv_intersection_agreement/svim/qual_15/", SV_INTERSECTION_HEADER)
write_results(pbsv_sv_intersection_means, "4_results/full_depth/sv_intersection_agreement/pbsv/", SV_INTERSECTION_HEADER)

# Subsampled Analysis - Summarize structural variant variant agreement for all strains
for depth in SUBSAMPLE_DEPTHS:
	
	subdepth = "subsampled/" + depth 
	sniffles_subsampled_path = SNIFFLES_PATH.replace("DEPTH", subdepth)
	#sniffles_subsampled_path = sniffles_subsampled_path.replace("/", subdepth, 1)
	sniffles_sv_intersection_means = get_agreement_summary("sniffles",sniffles_subsampled_path, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_total_svs.csv")
	outdir = "4_results/" + subdepth + "/sv_intersection_agreement/sniffles/"
	write_results(sniffles_sv_intersection_means, outdir, SV_INTERSECTION_HEADER)
	
	pbsv_subsampled_path = PBSV_PATH.replace("DEPTH", subdepth)
	pbsv_sv_intersection_means = get_agreement_summary("pbsv",pbsv_subsampled_path, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_total_svs.csv")
	outdir = "4_results/" + subdepth + "/sv_intersection_agreement/pbsv/"
	write_results(pbsv_sv_intersection_means, outdir, SV_INTERSECTION_HEADER)
	
	svim_subsampled_path = SVIM_PATH_QUAL_15.replace("DEPTH", subdepth)
	svim_sv_intersection_means = get_agreement_summary("svim",svim_subsampled_path, SNIFFLES_SVIM_ALIGNERS , SVIM_SV_TYPES, 1, "summary_total_svs.csv")
	outdir = "4_results/" + subdepth + "/sv_intersection_agreement/svim/qual_15/"
	write_results(svim_sv_intersection_means, outdir, SV_INTERSECTION_HEADER)	

# Summarize breakpoint agreement for all strains
sniffles_path = SNIFFLES_PATH.replace("DEPTH","full_depth")
sniffles_breakpoint_means = get_agreement_summary("sniffles",sniffles_path, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_breakpoints.csv")
svim_path = SVIM_PATH_QUAL_15.replace("DEPTH","full_depth")
svim_qual15_breakpoint_means = get_agreement_summary("svim", svim_path, SNIFFLES_SVIM_ALIGNERS, SVIM_SV_TYPES, 1, "summary_breakpoints.csv")
pbsv_path = PBSV_PATH.replace("DEPTH","full_depth")
pbsv_breakpoint_means = get_agreement_summary("pbsv",pbsv_path, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_breakpoints.csv")

write_results(sniffles_breakpoint_means, "4_results/full_depth/breakpoint_agreement/sniffles/", BREAKPOINT_HEADER)
write_results(svim_qual15_breakpoint_means, "4_results/full_depth/breakpoint_agreement/svim/qual_15/", BREAKPOINT_HEADER)
write_results(pbsv_breakpoint_means, "4_results/full_depth/breakpoint_agreement/pbsv/", BREAKPOINT_HEADER)

# Subsampled Analysis - Summarize breakpoint agreement for all strains
for depth in SUBSAMPLE_DEPTHS:
	subdepth = "subsampled/" + depth 
	sniffles_subsampled_path = SNIFFLES_PATH.replace("DEPTH",subdepth)
	#sniffles_subsampled_path = sniffles_subsampled_path.replace("/", subdepth, 1)
	sniffles_breakpoint_means = get_agreement_summary("sniffles",sniffles_subsampled_path, SNIFFLES_SVIM_ALIGNERS, SNIFFLES_SV_TYPES, 1, "summary_breakpoints.csv")
	outdir = "4_results/" + subdepth + "/breakpoint_agreement/sniffles/"
	write_results(sniffles_breakpoint_means, outdir, BREAKPOINT_HEADER)
	
	pbsv_subsampled_path = PBSV_PATH.replace("DEPTH",subdepth)
	pbsv_breakpoint_means = get_agreement_summary("pbsv",pbsv_subsampled_path, PBSV_ALIGNERS, PBSV_SV_TYPES, 1, "summary_breakpoints.csv")
	outdir = "4_results/" + subdepth + "/breakpoint_agreement/pbsv/"
	write_results(pbsv_breakpoint_means, outdir, BREAKPOINT_HEADER)
	
	svim_subsampled_path = SVIM_PATH_QUAL_15.replace("DEPTH",subdepth)
	svim_breakpoint_means = get_agreement_summary("svim",svim_subsampled_path, SNIFFLES_SVIM_ALIGNERS , SVIM_SV_TYPES, 1, "summary_breakpoints.csv")
	outdir = "4_results/" + subdepth + "/breakpoint_agreement/svim/qual_15/"
	
	write_results(svim_breakpoint_means, outdir, BREAKPOINT_HEADER)	

