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
import argparse
import sys
import csv

parser = argparse.ArgumentParser()
parser.add_argument("caller", help="Caller that generated the predicted SVs")
parser.add_argument("caller_csv_path", help="Path containing the caller csvs for each strain") # 3_variant_calls/DEPTH/ALIGNER/CALLER/STRAIN/shuffled/summary/
parser.add_argument("outputpath", help="Directory to write results to")
parser.add_argument('--svim_qual', default='15', help="SVIM qual value")
args = parser.parse_args()

ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "N2", "DL238","ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]
SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["pbmm2"]
PBSV_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV","ALL_SVS"]
SNIFFLES_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV","ALL_SVS"]
SVIM_SV_TYPES = ["BND", "DEL", "DUP:INT","DUP:TANDEM", "INS","INV","ALL_SVS"]

SV_INTERSECTION_HEADER = "SV TYPE,Intersection,Original only,Shuffled only,Unique,Unique proportion"
OVERLAP_SUMMARY_HEADER = "SV TYPE,Total calls,Intersection,Non-intersecting,Non-intersecting proportion"
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

		if line_count > 1:
			calls_total = 0
			intersection_total = 0
			non_intersecting_total = 0
			total_filters = []
			for line in sv_call_lines:
				total_filters.append(line[0])
				calls_total += int(line[1])
				intersection_total += int(line[2])
				non_intersecting_total += int(line[3])

			non_intersecting_proportion = non_intersecting_total / calls_total
			total_filters_joined = ";".join(total_filters)
			new_line = [strain, str(total_filters_joined), calls_total, intersection_total,non_intersecting_total,non_intersecting_proportion]
			sv_type_dict_totals[x].append(new_line)

		else:
			sv_agreement_line = list(chain.from_iterable(sv_call_lines))
			sv_agreement_line.insert(0,strain)
			sv_type_dict_totals[x].append(sv_agreement_line)

	return(sv_type_dict_totals)

# Function calculates the agreement statistics for each caller/aligner combo.
# Parameters: caller = caller name ,caller_csv_path = path to caller csv files, aligners = aligners used to generate alignments for the caller, caller_sv_types = sv types predicted by the caller
def get_agreement_summary(caller,caller_csv_path, aligners, caller_sv_types, inputfile, write_stdev):
	# caller_aligner_dict_bp: dictionary to store the sv overlap or breakpoint agreement data imported from csv files
	# caller_aligner_dict_bp: keys = caller/aligner combo and the svtype
	# caller_aligner_dict_bp: values = filters, intersection, same breakpoints, different breakpoints
	caller_aligner_dict_bp = defaultdict(lambda: defaultdict(list))

	# Loop is used to populate caller_aligner_dict_bp
	for aligner in aligners: # analyze results for each different aligner used for the caller
		caller_aligner_combo = caller + "-" + aligner
		for strain in ALL_STRAINS:
			caller_csv_file = caller_csv_path.replace("ALIGNER",aligner).replace("STRAIN",strain) + inputfile
			if os.path.isfile(caller_csv_file):
				# sv_agreement_dict stores csv file information
				# sv_agreement_dict: keys = Variant type, values = filters, intersection, same breakpoints or original only, different breakpoints or shuffled only, discordant proportion or total unique, unique proportion (only if summary_overlap_total_svs.csv was the csv file)
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
		if "overlap_comparison" in inputfile:
			for svtype in caller_aligner_dict_bp[caller_aligner]:
				total_calls_values = []
				intersection_values = []
				non_intersecting_values = []
				non_intersecting_proportion_values = []
				filters = set()

				for stat_line in caller_aligner_dict_bp[caller_aligner][svtype]:
					filters.add(stat_line[1])
					total_calls_values.append(int(stat_line[2]))
					intersection_values.append(int(stat_line[3]))
					non_intersecting_values.append(int(stat_line[4]))
					non_intersecting_proportion_values.append(float(stat_line[5]))

				mean_total_calls = statistics.mean(total_calls_values)
				mean_intersection = statistics.mean(intersection_values)
				mean_non_intersecting = statistics.mean(non_intersecting_values)
				mean_non_intersecting_proportion = statistics.mean(non_intersecting_proportion_values)
				stdev_total_calls = statistics.stdev(total_calls_values)
				stdev_intersection = statistics.stdev(intersection_values)
				stdev_non_intersecting = statistics.stdev(non_intersecting_values)
				stdev_non_intersecting_proportion = statistics.stdev(non_intersecting_proportion_values)
				if write_stdev:
					summary_line =f'{mean_total_calls:.0f}±{stdev_total_calls:.0f},{mean_intersection:.0f}±{stdev_intersection:.0f},{mean_non_intersecting:.0f}±{stdev_non_intersecting:.0f},{mean_non_intersecting_proportion:.3f}±{stdev_non_intersecting_proportion:.3f}'
				else:
					summary_line =f'{mean_total_calls:.0f},{mean_intersection:.0f},{mean_non_intersecting:.0f},{mean_non_intersecting_proportion:.3f}'
				caller_aligner_dict_bp_means[caller_aligner][svtype]= summary_line

		else:
			print("Wrong file name: " + inputfile)

	return(caller_aligner_dict_bp_means)


# Write results to disk
def write_results(result_dict, outputpath, header, filesuffix):
	if not os.path.exists(outputpath):
		os.makedirs(outputpath)
	for caller_aligner in result_dict.keys():
		outfile = outputpath + caller_aligner + filesuffix + ".csv"
		print("Writing to: " + outfile)
		csv_list = [header + "\n"] # List to store csv lines to be written.
		all_svs_line = None

		for svtype in sorted(result_dict[caller_aligner]):
			csv_line = svtype + "," + result_dict[caller_aligner][svtype] + "\n"

			if "ALL_SVS" in csv_line:
				all_svs_line = csv_line
			else:
				csv_list.append(csv_line)
		csv_list.append(all_svs_line)

		with open(outfile, 'w', newline='\n') as f:
			f.writelines(csv_list)


# Summarize structural variant overlap and breakpoint agreement for all strains
if args.caller == "pbsv":
	caller = "pbsv"
	aligners = PBSV_ALIGNERS
	caller_sv_types = PBSV_SV_TYPES

elif args.caller == "sniffles":
	caller = "sniffles"
	aligners = SNIFFLES_SVIM_ALIGNERS
	caller_sv_types = SNIFFLES_SV_TYPES

elif args.caller == "svim":
	caller = "svim"
	aligners = SNIFFLES_SVIM_ALIGNERS
	caller_sv_types = SVIM_SV_TYPES
	sv_intersection_outputdir = args.outputpath + "/sv_intersection_agreement/" + caller + "/qual_" + args.svim_qual + "/"

# Get means without standard deviations
sv_intersection_means = get_agreement_summary(caller,args.caller_csv_path, aligners, caller_sv_types , "overlap_comparison_all_svs.csv", False)
sv_intersection_means_coords = get_agreement_summary(caller,args.caller_csv_path, aligners, caller_sv_types , "overlap_comparison_coords_all_svs.csv", False)
sv_intersection_means_relaxed = get_agreement_summary(caller,args.caller_csv_path, aligners, caller_sv_types , "overlap_comparison_relaxed_all_svs.csv", False)

if args.caller != "svim":
	sv_intersection_outputdir = args.outputpath + "/sv_intersection_agreement/" + caller + "/"

write_results(sv_intersection_means, sv_intersection_outputdir, OVERLAP_SUMMARY_HEADER, "_agreement_summary_total")
write_results(sv_intersection_means_coords, sv_intersection_outputdir, OVERLAP_SUMMARY_HEADER, "_agreement_summary_coords_total")
write_results(sv_intersection_means_relaxed, sv_intersection_outputdir, OVERLAP_SUMMARY_HEADER, "_agreement_summary_relaxed_total")

# Get means and standard deviations
sv_intersection_means = get_agreement_summary(caller,args.caller_csv_path, aligners, caller_sv_types , "overlap_comparison_all_svs.csv", True)
sv_intersection_means_coords = get_agreement_summary(caller,args.caller_csv_path, aligners, caller_sv_types , "overlap_comparison_coords_all_svs.csv", True)
sv_intersection_means_relaxed = get_agreement_summary(caller,args.caller_csv_path, aligners, caller_sv_types , "overlap_comparison_relaxed_all_svs.csv", True)

if args.caller != "svim":
	sv_intersection_outputdir = args.outputpath + "/sv_intersection_agreement/" + caller + "/with_stdev/"
else:
	sv_intersection_outputdir = args.outputpath + "/sv_intersection_agreement/" + caller + "/qual_" + args.svim_qual + "/with_stdev/"

write_results(sv_intersection_means, sv_intersection_outputdir, OVERLAP_SUMMARY_HEADER, "_agreement_summary_total")
write_results(sv_intersection_means_coords, sv_intersection_outputdir, OVERLAP_SUMMARY_HEADER, "_agreement_summary_coords_total")
write_results(sv_intersection_means_relaxed, sv_intersection_outputdir, OVERLAP_SUMMARY_HEADER, "_agreement_summary_relaxed_total")
