import argparse
import os
import collections
from collections import defaultdict
import pandas
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("original_vcf", help="Original vcf")
parser.add_argument("shuffled_vcf", help="Shuffled vcf")
parser.add_argument("output_dir", help="Output directory")
parser.add_argument("sv_caller", help="variant caller used to create vcf file. Use pbsv, sniffles, or svim")
parser.add_argument("--minsize", default=0, type=int, help="minimum variant size in bp")
parser.add_argument("--min_qual_svim", default=0, type=int, help="minimum qual value for svim")
parser.add_argument("--reciprocal", type=str, help="minimum reciprocal overlap for calls in the original and shuffled sets to be considered as Overlapping")

args = parser.parse_args()

# Types of variants to look for
PBSV_SV_TYPES = ["INV","DEL","DUP", "INS", "BND"]
SNIFFLES_SV_TYPES = ["INV", "DEL", "DUP", "INVDUP", "INV/INVDUP", "DEL/INV", "TRA", "INS", "BND", "DUP/INS"] # Note: DUP/INS may cause problems, as some have END coords that start before the start. DUP/INS is rare, so they can probably be ignored altogether.
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

# Parse VCF file for given caller
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

# Create dictionary for SVs in all VCF files
def create_dict_all_vcfs(all_dicts):
	total_sv_types = set()
	sv_dict_all_vcfs = defaultdict(set) # Store all variant coordinates in dictionary using filters and sv types as keys
	for vcf_dict in all_dicts:
		for sv_type in vcf_dict.keys():
			total_sv_types.add(sv_type)
			for sv_call in vcf_dict[sv_type]:
				sv_dict_all_vcfs[sv_type].add(sv_call)
	return(sv_dict_all_vcfs, total_sv_types)

# Summarize VCF file and store in log
def summarize_vcf(total_svs_1, vcf_file):
	summary_log = []
	summary_log.append(("Results for : " + vcf_file))

	for sv_type in total_svs_1.keys():
		sv_call_count = len(total_svs_1[sv_type])
		summary_log.append(("\t" + sv_type + ": " + str(sv_call_count)))

	return(summary_log)

def get_vcf_intersection(all_dicts):
	vcf_file_log = []
	sv_dict_1 = all_dicts[0]
	sv_dict_1_sv_types = set(sv_dict_1.keys())
	sv_dict_2 = all_dicts[1]
	sv_dict_2_sv_types = set(sv_dict_2.keys())

	sv_type_intersection = set.intersection(sv_dict_1_sv_types, sv_dict_2_sv_types)
	sv_type_union = set.union(sv_dict_1_sv_types, sv_dict_2_sv_types)
	sv_type_difference = set.difference(sv_type_union, sv_type_intersection)

	if len(sv_type_difference) > 0:
		print("Warning: Not all sv types found in all VCF files")
		print(sv_type_difference)
	else:
		sv_intersection_stats_original_shuffled_df = pandas.DataFrame(columns=['Union', 'Intersection', 'Symmetric difference'])
		for sv_type in sv_type_intersection:
			sv_dict_1_calls = sv_dict_1[sv_type]
			sv_dict_2_calls = sv_dict_2[sv_type]

			# Get intersection, union, and difference
			sv_call_intersection = set.intersection(sv_dict_1_calls, sv_dict_2_calls)
			sv_call_union = set.union(sv_dict_1_calls,sv_dict_2_calls)
			sv_call_differences = set.difference(sv_call_union, sv_call_intersection) # Symmetric difference
			vcf_file_log.append("SV Type: " + sv_type)
			vcf_file_log.append("\tIntersection: " + str(len(sv_call_intersection)))
			vcf_file_log.append("\tSymmetric difference: " + str(len(sv_call_differences)))

			sv_intersection_stats_original_shuffled = {'Union': len(sv_call_union), 'Intersection': len(sv_call_intersection), 'Symmetric difference': len(sv_call_differences)}
			sv_intersection_stats_original_shuffled_series = pandas.Series(sv_intersection_stats_original_shuffled, name = sv_type)
			sv_intersection_stats_original_shuffled_df_new_row = pandas.DataFrame([sv_intersection_stats_original_shuffled_series], columns=sv_intersection_stats_original_shuffled_df.columns)
			sv_intersection_stats_original_shuffled_df = pandas.concat([sv_intersection_stats_original_shuffled_df, sv_intersection_stats_original_shuffled_df_new_row])

	sv_intersection_stats_original_shuffled_df['Union'] = pandas.to_numeric(sv_intersection_stats_original_shuffled_df['Union'])
	sv_intersection_stats_original_shuffled_df['Intersection'] = pandas.to_numeric(sv_intersection_stats_original_shuffled_df['Intersection'])
	sv_intersection_stats_original_shuffled_df['Symmetric difference'] = pandas.to_numeric(sv_intersection_stats_original_shuffled_df['Symmetric difference'])
	sv_intersection_stats_original_shuffled_df.loc['Total']= sv_intersection_stats_original_shuffled_df.sum(numeric_only=True, axis=0)

	return(vcf_file_log, sv_intersection_stats_original_shuffled_df)

from pathlib import Path


def write_lines(my_list, outfile):
	with open(outfile, 'w', newline='\n') as f:
		for line in my_list:
			f.write(f"{line}\n")

def get_intersection_proportions(my_df):
	jaccard_df = my_df[['SV Type']].copy()
	jaccard_df['Jaccard Index'] = pandas.to_numeric(my_df['Intersection'] / my_df['Union'])
	jaccard_df['Jaccard Distance'] = pandas.to_numeric(my_df['Symmetric Difference'] / my_df['Union'])
	jaccard_df['Symmetric Difference'] = pandas.to_numeric(my_df['Symmetric Difference'])

	return(jaccard_df)

vcf_file_1 = parse_vcf_file(args.original_vcf)
total_svs_1 = vcf_file_1[0]
total_svs_1_no_meta = vcf_file_1[1]

vcf_file_2 = parse_vcf_file(args.shuffled_vcf)
total_svs_2 = vcf_file_2[0]
total_svs_2_no_meta = vcf_file_2[1]

sv_dict_all_vcfs = create_dict_all_vcfs([total_svs_1, total_svs_2])
summary_log = summarize_vcf(total_svs_1,args.original_vcf)
summary_log.extend(summarize_vcf(total_svs_2,args.shuffled_vcf))

Path(args.output_dir).mkdir(parents=True, exist_ok=True)
logfile = args.output_dir + "/vcf_summary.log"
write_lines(summary_log, logfile)

summarize_vcf(sv_dict_all_vcfs[0], "all vcfs")

# Calculate intersection statistics for SV calls including all metadata
vcf_file_intersections = get_vcf_intersection([total_svs_1, total_svs_2])
vcf_file_intersections_log = vcf_file_intersections[0]
vcf_file_intersections_results = vcf_file_intersections[1]
logfile = args.output_dir + "/sv_intersection_stats.log"

write_lines(vcf_file_intersections_log, logfile)

vcf_file_intersections_results = vcf_file_intersections_results.reset_index()
vcf_file_intersections_results.rename(columns={'index':'SV Type'}, inplace=True)
vcf_file_intersections_results.sort_values(by=['SV Type'], inplace=True)

vcf_file_intersections_results_prop = get_intersection_proportions(vcf_file_intersections_results)
results_file = args.output_dir + "/results.csv"
vcf_file_intersections_results.to_csv(results_file, index = False)
prop_results_file = args.output_dir + "/results_proportions.csv"

vcf_file_intersections_results_prop.to_csv(prop_results_file, index = False, float_format='%f')

# Calculate intersection statistics for SV calls excluding sequences, INFO, and FORMAT fields (coordinates only)
vcf_file_intersections_no_meta = get_vcf_intersection([total_svs_1_no_meta, total_svs_2_no_meta])
vcf_file_intersections_no_meta_log = vcf_file_intersections_no_meta[0]
vcf_file_intersections_no_meta_results = vcf_file_intersections_no_meta[1]
logfile = args.output_dir + "/sv_intersection_stats_no_meta.log"

write_lines(vcf_file_intersections_no_meta_log, logfile)

vcf_file_intersections_no_meta_results = vcf_file_intersections_no_meta_results.reset_index()
vcf_file_intersections_no_meta_results.rename(columns={'index':'SV Type'}, inplace=True)
vcf_file_intersections_no_meta_results.sort_values(by=['SV Type'], inplace=True)

vcf_file_intersections_no_meta_results_prop = get_intersection_proportions(vcf_file_intersections_no_meta_results)
results_file = args.output_dir + "/results_no_meta.csv"
vcf_file_intersections_no_meta_results.to_csv(results_file, index = False)
prop_results_file = args.output_dir + "/results_proportions_no_meta.csv"
vcf_file_intersections_no_meta_results_prop.to_csv(prop_results_file, index = False, float_format='%f')
