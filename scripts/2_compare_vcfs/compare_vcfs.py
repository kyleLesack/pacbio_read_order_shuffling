import argparse
import os
import collections
from collections import defaultdict
import pandas
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("input_vcf_1", help="1st vcf file to compare")
parser.add_argument("input_vcf_2", help="2nd vcf file to compare")
parser.add_argument("input_vcf_3", help="3rd vcf file to compare")
parser.add_argument("input_vcf_4", help="4th vcf file to compare")
parser.add_argument("input_vcf_5", help="5th vcf file to compare")
parser.add_argument("input_vcf_6", help="6th vcf file to compare")
parser.add_argument("output_dir", help="Output directory")

parser.add_argument("sv_caller", help="variant caller used to create vcf file. Use pbsv, sniffles, or svim")
parser.add_argument("--minsize", default=0, type=int, help="minimum variant size in bp")
parser.add_argument("--min_qual_svim", default=0, type=int, help="minimum qual value for svim")
parser.add_argument("--reciprocal", type=str, help="minimum reciprocal overlap for calls in the original and shuffled sets to be considered as Overlapping")

args = parser.parse_args()

# Types of variants to look for
PBSV_SV_TYPES = ["INV","DEL","DUP", "INS", "BND"]
SNIFFLES_SV_TYPES = ["INV", "DEL", "DUP", "INVDUP", "INV/INVDUP", "DEL/INV", "TRA", "INS", "BND", "DUP/INS"] # Note: DUP/INS causes problems, as some have END coords that start before the start. DUP/INS is rare, so I can probably ignore them altogether.
SVIM_SV_TYPES = ["INV", "DEL","DUP:INT", "DUP:TANDEM","INS","DUP","BND"]

# Extract pbsv variants
def parse_pbsv(pbsv_variants):
	sv_dict_total_vcf_line = defaultdict(set) # Store all vcf lines in dictionary using sv types as keys
	sv_dict_total_no_meta = defaultdict(set) # Store all vcf sv calls in dictionary sv types as keys. Only store CHR, START, END, SEQ, FILTER
	excluded_due_to_type = set()
	excluded_due_to_size = set()

	# Go through SVs and extract desired information
	for line in pbsv_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			variant_seq = line_split[4] # SV sequence for most SV types. BND describes where the other coord mapped to (e.g., translocation). See https://github.com/eldariont/pbsv/issues/28
			line = line.replace(variant_name, "pbsv_sv")
			filter_line = line_split[6] # Describes if variant call passed. Most calls are designated as PASSED, so I probably can ignore the rest, which are probably low quality anyways.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			end_coord = None
			variant_type = None
			variant_size = None
			variant_size_vcf = None # SV size stored in VCF file line

			# Get the end coordinate, variant type, and size
			for x in info_line_split:
				if "END=" in x:
					end_coord = x.split("=")[1]
					variant_size = int(end_coord) - int(start_coord)
				elif "SVTYPE=" in x:
					variant_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = int(x.split("=")[1])
			if variant_size_vcf is not None:
				variant_size = variant_size_vcf

			if variant_type is not None:
				if variant_type not in PBSV_SV_TYPES:
					excluded_due_to_type.add(line)
				else:
					if variant_type != "BND":
						if variant_size is None:
							print(line)
							input("missing sv size")

						else:
							abs_variant_size= abs(int(variant_size))

						if variant_type in PBSV_SV_TYPES and abs_variant_size >= args.minsize:
							if variant_type == "INS":
								sv_dict_total_vcf_line[variant_type].add(line)
								sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + variant_seq
								sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
							else:
								sv_dict_total_vcf_line[variant_type].add(line)
								sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_seq
								sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
						else:
							if variant_type not in PBSV_SV_TYPES:
								print("Excluding: " + line)
								print("Invalid variant type: " + variant_type)
							else:
								excluded_due_to_size.add(line)
					else:
						sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + variant_seq
						sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
						sv_dict_total_vcf_line[variant_type].add(line)
			else:
				print("Variant type is None")
				print(line)

	return (sv_dict_total_vcf_line, sv_dict_total_no_meta, excluded_due_to_type, excluded_due_to_size)

# Extract sniffles variants
def parse_sniffles(sniffles_variants):
	sv_dict_total_vcf_line = defaultdict(set) # Store all vcf lines in dictionary using sv types as keys
	sv_dict_total_no_meta = defaultdict(set) # Store all vcf sv calls in dictionary sv types as keys. Only store CHR, START, END, SEQ, FILTER
	excluded_due_to_type = set()
	excluded_due_to_size = set()

	# Go through SVs and extract desired information
	for line in sniffles_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			variant_seq = line_split[4] # SV sequence for most SV types. BND describes where the other coord mapped to (e.g., translocation). See https://github.com/eldariont/sniffles/issues/28
			line = line.replace(variant_name, "sniffles_sv")
			filter_line = line_split[6] # Describes if variant call passed. Most calls are designated as PASSED, so I probably can ignore the rest, which are probably low quality anyways.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			end_coord = None
			variant_type = None
			variant_size = None
			variant_size_vcf = None # SV size stored in VCF file line

			# Get the end coordinate, variant type, and size
			for x in info_line_split:
				if "END=" in x:
					end_coord = x.split("=")[1]
					variant_size = int(end_coord) - int(start_coord)
				elif "SVTYPE=" in x:
					variant_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = int(x.split("=")[1])
			if variant_size_vcf is not None:
				variant_size = variant_size_vcf

			if variant_type is not None:
				if variant_type not in SNIFFLES_SV_TYPES:
					excluded_due_to_type.add(line)
				else:
					if variant_type != "BND":
						if variant_size is None:
							print(line)
							input("missing sv size")

						else:
							abs_variant_size= abs(int(variant_size))

						if variant_type in SNIFFLES_SV_TYPES and abs_variant_size >= args.minsize:
							if variant_type == "INS":
								sv_dict_total_vcf_line[variant_type].add(line)
								sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + variant_seq
								sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
							else:
								sv_dict_total_vcf_line[variant_type].add(line)
								sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_seq
								sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
						else:
							if variant_type not in SNIFFLES_SV_TYPES:
								print("Excluding: " + line)
								print("Invalid variant type: " + variant_type)
							else:
								excluded_due_to_size.add(line)
					else:
						sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + variant_seq
						sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
						sv_dict_total_vcf_line[variant_type].add(line)
			else:
				print("Variant type is None")
				print(line)

	return (sv_dict_total_vcf_line, sv_dict_total_no_meta, excluded_due_to_type, excluded_due_to_size)

# Extract svim variants
def parse_svim(svim_variants):
	sv_dict_total_vcf_line = defaultdict(set) # Store all vcf lines in dictionary using sv types as keys
	sv_dict_total_no_meta = defaultdict(set) # Store all vcf sv calls in dictionary sv types as keys. Only store CHR, START, END, SEQ, FILTER
	excluded_due_to_type = set()
	excluded_due_to_size = set()

	# Go through SVs and extract desired information
	for line in svim_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			variant_seq = line_split[4] # SV sequence for most SV types. BND describes where the other coord mapped to (e.g., translocation). See https://github.com/eldariont/svim/issues/28
			line = line.replace(variant_name, "svim_sv")
			variant_qual = int(line_split[5]) # svim qual score
			filter_line = line_split[6] # Describes if variant call passed. Most calls are designated as PASSED, so I probably can ignore the rest, which are probably low quality anyways.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			end_coord = None
			variant_type = None
			variant_size = None
			variant_size_vcf = None # SV size stored in VCF file line

			# Get the end coordinate, variant type, and size
			for x in info_line_split:
				if "END=" in x:
					end_coord = x.split("=")[1]
					variant_size = int(end_coord) - int(start_coord)
				elif "SVTYPE=" in x:
					variant_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = int(x.split("=")[1])
			if variant_size_vcf is not None:
				variant_size = variant_size_vcf

			if variant_type is not None:
				if variant_type not in SVIM_SV_TYPES:
					excluded_due_to_type.add(line)
				else:
					if variant_type != "BND":
						if variant_size is None:
							print(line)
							input("missing sv size")

						else:
							abs_variant_size= abs(int(variant_size))

						if variant_type in SVIM_SV_TYPES and variant_qual >= args.min_qual_svim and abs_variant_size >= args.minsize:
							if variant_type == "INS":
								sv_dict_total_vcf_line[variant_type].add(line)
								sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + variant_seq
								sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
							else:
								sv_dict_total_vcf_line[variant_type].add(line)
								sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_seq
								sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
						else:
							if variant_type not in SVIM_SV_TYPES:
								print("Excluding: " + line)
								print("Invalid variant type: " + variant_type)
							else:
								excluded_due_to_size.add(line)
					else:
						if variant_qual >= args.min_qual_svim:
							sv_call_no_meta = filter_line + "\t" + chromosome + "\t" + start_coord + "\t" + "\t" + variant_seq
							sv_dict_total_no_meta[variant_type].add(sv_call_no_meta)
							end_coord = line_split[4]
							sv_dict_total_vcf_line[variant_type].add(line)
			else:
				print("Variant type is None")
				print(line)

	return (sv_dict_total_vcf_line, sv_dict_total_no_meta, excluded_due_to_type, excluded_due_to_size)

# Parse vcf file for given caller
def parse_vcf_file(vcf_file):
	with open(vcf_file) as f:
		variants_1 = f.readlines()
		if args.sv_caller.lower() == "pbsv":
			vcf_file_1 = parse_pbsv(variants_1)
		elif args.sv_caller.lower() == "sniffles":
			vcf_file_1 = parse_sniffles(variants_1)
		elif args.sv_caller.lower() == "svim":
			vcf_file_1 = parse_svim(variants_1)

	return vcf_file_1

# Create a dictionary with svs from all vcf files
def create_dict_all_vcfs(all_dicts):
	total_sv_types = set()
	sv_dict_all_vcfs = defaultdict(set) # Store all variant coordinates in dictionary using filters and sv types as keys
	for vcf_dict in all_dicts:
		for sv_type in vcf_dict.keys():
			total_sv_types.add(sv_type)
			for sv_call in vcf_dict[sv_type]:
				sv_dict_all_vcfs[sv_type].add(sv_call)
	return(sv_dict_all_vcfs, total_sv_types)

# Summarize the vcf file and save in log file
def summarize_vcf(total_svs_1, vcf_file):
	summary_log = []
	summary_log.append(("Results for : " + vcf_file))

	for sv_type in total_svs_1.keys():
		sv_call_count = len(total_svs_1[sv_type])
		summary_log.append(("\t" + sv_type + ": " + str(sv_call_count)))

	return(summary_log)

# Get agreement between different vcfs
def get_vcf_intersection(all_dicts):
	vcf_file_log = []
	sv_dict_1 = all_dicts[0]
	sv_dict_1_sv_types = set(sv_dict_1.keys())
	sv_dict_2 = all_dicts[1]
	sv_dict_2_sv_types = set(sv_dict_2.keys())
	sv_dict_3 = all_dicts[2]
	sv_dict_3_sv_types = set(sv_dict_3.keys())
	sv_dict_4 = all_dicts[3]
	sv_dict_4_sv_types = set(sv_dict_4.keys())
	sv_dict_5 = all_dicts[4]
	sv_dict_5_sv_types = set(sv_dict_5.keys())
	sv_dict_6 = all_dicts[5]
	sv_dict_6_sv_types = set(sv_dict_6.keys())

	# Get intersection, union, and differences
	sv_type_intersection = set.intersection(sv_dict_1_sv_types, sv_dict_2_sv_types, sv_dict_3_sv_types, sv_dict_4_sv_types, sv_dict_5_sv_types, sv_dict_6_sv_types)
	sv_type_union = set.union(sv_dict_1_sv_types, sv_dict_2_sv_types, sv_dict_3_sv_types, sv_dict_4_sv_types, sv_dict_5_sv_types, sv_dict_6_sv_types)
	sv_type_difference = set.difference(sv_type_union, sv_type_intersection)

	if len(sv_type_difference) > 0:
		print("Warning: Not all sv types found in all VCF files")
		print(sv_type_difference)
	else:
		num_callers_with_variant_df = pandas.DataFrame(columns=['One', 'Two', 'Three', 'Four', 'Five', 'Six'])
		for sv_type in sv_type_intersection:
			sv_dict_1_calls = sv_dict_1[sv_type]
			sv_dict_2_calls = sv_dict_2[sv_type]
			sv_dict_3_calls = sv_dict_3[sv_type]
			sv_dict_4_calls = sv_dict_4[sv_type]
			sv_dict_5_calls = sv_dict_5[sv_type]

			sv_call_intersection = set.intersection(sv_dict_1_calls, sv_dict_2_calls, sv_dict_3_calls, sv_dict_4_calls, sv_dict_5_calls)
			sv_call_union = set.union(sv_dict_1_calls,sv_dict_2_calls,sv_dict_3_calls,sv_dict_4_calls,sv_dict_5_calls)
			sv_call_differences = set.difference(sv_call_union, sv_call_intersection)
			vcf_file_log.append("SV Type: " + sv_type)
			vcf_file_log.append("\tIntersection: " + str(len(sv_call_intersection)))
			vcf_file_log.append("\tDifferences: " + str(len(sv_call_differences)))

			callers_with_sv_dict = {'One': 0, 'Two': 0, 'Three': 0, 'Four': 0, 'Five':0, 'Six': 0}

			for sv_call in sv_call_union:
				sv_callers_with_sv_call = 0
				if sv_call in sv_dict_1[sv_type]:
					sv_callers_with_sv_call +=1
				if sv_call in sv_dict_2[sv_type]:
					sv_callers_with_sv_call +=1
				if sv_call in sv_dict_3[sv_type]:
					sv_callers_with_sv_call +=1
				if sv_call in sv_dict_4[sv_type]:
					sv_callers_with_sv_call +=1
				if sv_call in sv_dict_5[sv_type]:
					sv_callers_with_sv_call +=1
				if sv_call in sv_dict_6[sv_type]:
					sv_callers_with_sv_call +=1

				# Update dictionary to reflect number of callers with the same sv call
				if sv_callers_with_sv_call == 1:
					callers_with_sv_dict['One'] +=1
				elif sv_callers_with_sv_call == 2:
					callers_with_sv_dict['Two'] +=1
				elif sv_callers_with_sv_call == 3:
					callers_with_sv_dict['Three'] +=1
				elif sv_callers_with_sv_call == 4:
					callers_with_sv_dict['Four'] +=1
				elif sv_callers_with_sv_call == 5:
					callers_with_sv_dict['Five'] +=1
				elif sv_callers_with_sv_call == 6:
					callers_with_sv_dict['Six'] +=1
				else:
					input("Invalid value for sv_callers_with_sv_call")

			callers_with_sv_series = pandas.Series(callers_with_sv_dict, name = sv_type)

			num_callers_with_variant_df_new_row = pandas.DataFrame([callers_with_sv_series], columns=num_callers_with_variant_df.columns)
			num_callers_with_variant_df = pandas.concat([num_callers_with_variant_df, num_callers_with_variant_df_new_row])

	return(vcf_file_log, num_callers_with_variant_df)

from pathlib import Path


def write_lines(my_list, outfile):
	with open(outfile, 'w', newline='\n') as f:
		for line in my_list:
			f.write(f"{line}\n")

def get_intersection_proportions(my_df):
	my_df['Total'] =  my_df[['One', 'Two', 'Three', 'Four', 'Five', 'Six']].sum(axis=1)
	new_df = my_df[['SV Type']].copy()
	new_df['One'] = pandas.to_numeric(my_df['One'] / my_df['Total'])
	new_df['Two'] = pandas.to_numeric(my_df['Two'] / my_df['Total'])
	new_df['Three'] = pandas.to_numeric(my_df['Three'] / my_df['Total'])
	new_df['Four'] = pandas.to_numeric(my_df['Four'] / my_df['Total'])
	new_df['Five'] = pandas.to_numeric(my_df['Five'] / my_df['Total'])
	new_df['Six'] = pandas.to_numeric(my_df['Six'] / my_df['Total'])

	return(new_df)

vcf_file_1 = parse_vcf_file(args.input_vcf_1)
total_svs_1 = vcf_file_1[0]
total_svs_1_no_meta = vcf_file_1[1]

vcf_file_2 = parse_vcf_file(args.input_vcf_2)
total_svs_2 = vcf_file_2[0]
total_svs_2_no_meta = vcf_file_2[1]

vcf_file_3 = parse_vcf_file(args.input_vcf_3)
total_svs_3 = vcf_file_3[0]
total_svs_3_no_meta = vcf_file_3[1]

vcf_file_4 = parse_vcf_file(args.input_vcf_4)
total_svs_4 = vcf_file_4[0]
total_svs_4_no_meta = vcf_file_4[1]

vcf_file_5 = parse_vcf_file(args.input_vcf_5)
total_svs_5 = vcf_file_5[0]
total_svs_5_no_meta = vcf_file_5[1]

vcf_file_6 = parse_vcf_file(args.input_vcf_6)
total_svs_6 = vcf_file_6[0]
total_svs_6_no_meta = vcf_file_6[1]

sv_dict_all_vcfs = create_dict_all_vcfs([total_svs_1, total_svs_2, total_svs_3, total_svs_4, total_svs_5, total_svs_6])

summary_log = summarize_vcf(total_svs_1,args.input_vcf_1)
summary_log.extend(summarize_vcf(total_svs_2,args.input_vcf_2))
summary_log.extend(summarize_vcf(total_svs_3,args.input_vcf_3))
summary_log.extend(summarize_vcf(total_svs_4,args.input_vcf_4))
summary_log.extend(summarize_vcf(total_svs_5,args.input_vcf_5))
summary_log.extend(summarize_vcf(total_svs_6,args.input_vcf_6))

Path(args.output_dir).mkdir(parents=True, exist_ok=True)
logfile = args.output_dir + "/vcf_summary.log"
write_lines(summary_log, logfile)

summarize_vcf(sv_dict_all_vcfs[0], "all vcfs")

# Calculate intersection statistics for SV calls including all metadata
vcf_file_intersections = get_vcf_intersection([total_svs_1, total_svs_2, total_svs_3, total_svs_4, total_svs_5, total_svs_6])
vcf_file_intersections_log = vcf_file_intersections[0]
vcf_file_intersections_results = vcf_file_intersections[1]
logfile = args.output_dir + "/sv_intersection_summary.log"
write_lines(vcf_file_intersections_log, logfile)

vcf_file_intersections_results = vcf_file_intersections_results.reset_index()
vcf_file_intersections_results.rename(columns={'index':'SV Type'}, inplace=True)
vcf_file_intersections_results.sort_values(by=['SV Type'], inplace=True)

vcf_file_intersections_results_prop = get_intersection_proportions(vcf_file_intersections_results)
results_file = args.output_dir + "/results.csv"
vcf_file_intersections_results.to_csv(results_file, index = False)
prop_results_file = args.output_dir + "/results_proportions.csv"
vcf_file_intersections_results_prop.to_csv(prop_results_file, index = False, float_format='%.3f')

# Calculate intersection statistics for SV calls excluding sequences, INFO, and FORMAT fields (coordinates only)
vcf_file_intersections_no_meta = get_vcf_intersection([total_svs_1_no_meta, total_svs_2_no_meta, total_svs_3_no_meta, total_svs_4_no_meta, total_svs_5_no_meta, total_svs_6_no_meta])
vcf_file_intersections_no_meta_log = vcf_file_intersections_no_meta[0]
vcf_file_intersections_no_meta_results = vcf_file_intersections_no_meta[1]
logfile = args.output_dir + "/sv_intersection_summary_no_meta.log"
write_lines(vcf_file_intersections_no_meta_log, logfile)

vcf_file_intersections_no_meta_results = vcf_file_intersections_no_meta_results.reset_index()
vcf_file_intersections_no_meta_results.rename(columns={'index':'SV Type'}, inplace=True)
vcf_file_intersections_no_meta_results.sort_values(by=['SV Type'], inplace=True)

vcf_file_intersections_no_meta_results_prop = get_intersection_proportions(vcf_file_intersections_no_meta_results)
results_file = args.output_dir + "/results_no_meta.csv"
vcf_file_intersections_no_meta_results.to_csv(results_file, index = False)
prop_results_file = args.output_dir + "/results_proportions_no_meta.csv"
vcf_file_intersections_no_meta_results_prop.to_csv(prop_results_file, index = False, float_format='%.3f')
