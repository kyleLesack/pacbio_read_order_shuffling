import argparse
import os
import collections
from pybedtools import BedTool

parser = argparse.ArgumentParser()
parser.add_argument("shuffled_vcf", help="vcf file created from shuffled FASTQ")
parser.add_argument("original_vcf", help="vcf file created from original FASTQ") # The original FASTQ was used to create shuffled FASTQ files
parser.add_argument("sv_caller", help="variant caller used to create vcf file. Use pbsv, sniffles, or svim")
parser.add_argument("--minsize", default=0, type=int, help="minimum variant size in bp")
parser.add_argument("--min_qual_svim", default=0, type=int, help="minimum qual value for svim")
parser.add_argument("--reciprocal", type=str, help="minimum reciprocal overlap for calls in the original and shuffled sets to be considered as Overlapping")
args = parser.parse_args()

###TODO: Remove the BND code in the parse functions###

# Types of variants to look for
#PBSV_SV_TYPES = ["INV","DEL","DUP", "INS", "BND"]
PBSV_SV_TYPES = ["INV","DEL","DUP", "INS"]
#SNIFFLES_SV_TYPES = ["INV", "DEL", "DUP", "INVDUP", "INV/INVDUP", "DEL/INV", "TRA", "INS", "BND", "DUP/INS"] # Note: DUP/INS causes problems, as some have END coords that start before the start. DUP/INS is rare, so I can ignore them altogether.
SNIFFLES_SV_TYPES = ["INV", "DEL", "DUP", "INVDUP", "INV/INVDUP", "DEL/INV", "TRA", "INS", "DUP/INS"] # Note: DUP/INS causes problems, as some have END coords that start before the start. DUP/INS is rare, so I can ignore them altogether.
#SVIM_SV_TYPES = ["INV", "DEL","DUP:INT", "DUP:TANDEM","INS","DUP","BND"]
SVIM_SV_TYPES = ["INV", "DEL","DUP:INT", "DUP:TANDEM","INS","DUP"]

from collections import defaultdict

# Extract pbsv variants
def parse_pbsv(pbsv_variants):
	sv_dict_total = defaultdict(lambda: defaultdict(list)) # Store all variant coordinates in dictionary using filters and sv types as keys
	sv_dict_total_vcf_line = defaultdict(lambda: defaultdict(list)) # Store all vcf lines in dictionary using filters and sv types as keys
	#excluded_sv_types = set()
	excluded_due_to_size = []
	excluded_due_to_type = []

	for line in pbsv_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			line = line.replace(variant_name, "pbsv_sv")
			filter_line = line_split[6] # Describes if variant call passed.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			if info_line_split[0] == "IMPRECISE":
				filter_line = filter_line + "-IMPRECISE"
			end_coord = None
			variant_type = None
			variant_size = None

			for x in info_line_split: # Get the end coordinate, variant type, and size
				if "END=" in x:
					end_coord = x.split("=")[1]
					variant_size = int(end_coord) - int(start_coord)
				elif "SVTYPE=" in x:
					variant_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = x.split("=")[1]

			if variant_type is not None:
				if variant_type not in PBSV_SV_TYPES:
					excluded_due_to_type.append(line)
				else:
					if variant_type != "BND":
						if variant_type == "INS":
							abs_variant_size= abs(int(variant_size_vcf))
						else:
							variant_size_from_coords = int(end_coord) - int(start_coord)
							abs_variant_size= abs(int(variant_size_from_coords))

						if variant_type in PBSV_SV_TYPES and abs_variant_size >= args.minsize:
							if variant_type == "INS":
								bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
								sv_dict_total[filter_line][variant_type].append(bed_line)
								sv_dict_total_vcf_line[filter_line][variant_type].append(line)
							else:
								bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
								sv_dict_total[filter_line][variant_type].append(bed_line)
								sv_dict_total_vcf_line[filter_line][variant_type].append(line)

						else:
							if variant_type not in PBSV_SV_TYPES:
								print("Invalid variant type: " + variant_type)
							else:
								excluded_due_to_size.append(line)
					else:
						end_coord = line_split[4]
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						sv_dict_total[filter_line][variant_type].append(bed_line)
						sv_dict_total_vcf_line[filter_line][variant_type].append(line)
			else:
				print("Variant type is none")
				print(line)

	return (sv_dict_total_vcf_line, sv_dict_total,excluded_due_to_type, excluded_due_to_size)

# Extract variant calls from sniffles
def parse_sniffles(sniffles_variants):
	sv_dict_total = defaultdict(lambda: defaultdict(list)) # Store all variant coordinates in dictionary using filters and sv types as keys
	sv_dict_total_vcf_line = defaultdict(lambda: defaultdict(list)) # Store all vcf lines in dictionary using filters and sv types as keys
	excluded_sv_types = set()
	excluded_due_to_size = []
	excluded_due_to_type = []

	for line in sniffles_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			line = line.replace(variant_name, "sniffles_sv")
			filter_line = line_split[6] # Describes if variant call passed.
			info_line = line_split[7]
			info_line_split = info_line.split(";")
			variant_precision = info_line_split[0] # Describes if breakpoints are precise
			filter_line = filter_line + "-" + variant_precision # Combine filter with precision

			end_coord = None
			variant_type = None
			variant_size_vcf = None

			for x in info_line_split: # Get the end coordinate, variant type, and size
				if "END=" in x:
					end_coord = x.split("=")[1]
				elif "SVTYPE=" in x:
					variant_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = x.split("=")[1]

			if variant_type is not None:

				if variant_type not in SNIFFLES_SV_TYPES:
					excluded_due_to_type.append(line)
				else:

					variant_renamed = "sniffles_" + variant_type + "_" + variant_name
					support_line = line_split[9] # Line that provides read support for the reference allele and sv.
					sv_read_support = support_line.split(":")[2]

					if variant_type == "BND":
						end_coord = line_split[4]
						if end_coord is not None:
							if variant_type in SNIFFLES_SV_TYPES:
								bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
								sv_dict_total[filter_line][variant_type].append(bed_line)
								sv_dict_total_vcf_line[filter_line][variant_type].append(line)
							else:
								print("Excluding variant type: " + variant_type + " in vcf file: " + str(args.input_vcf))
						else:
							print("Excluding: " + line)
							if end_coord is None:
								print("End coordinate is None")

					elif variant_type is not None and end_coord is not None and variant_size_vcf is not None:
						if variant_type in SNIFFLES_SV_TYPES:
							abs_variant_size= abs(int(variant_size_vcf))

							if abs_variant_size >= args.minsize: # Select variants of a minimum size
								if variant_type == "INS":
									bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
									sv_dict_total[filter_line][variant_type].append(bed_line)
									sv_dict_total_vcf_line[filter_line][variant_type].append(line)
								else:
									bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
									sv_dict_total[filter_line][variant_type].append(bed_line)
									sv_dict_total_vcf_line[filter_line][variant_type].append(line)
							else:
								excluded_due_to_size.append(line)
						else:
							print("Excluding variant type: " + variant_type + " in vcf file: " + args.input_vcf)
					else:
						print("Excluding: " + line)
						if variant_type is None:
							print("Variant type is None")
						if end_coord is None:
							print("End coordinate is None")
						if variant_size_vcf is None:
							print("Variant size is None")

	return (sv_dict_total_vcf_line, sv_dict_total,excluded_due_to_type, excluded_due_to_size)

# Extract svim variants
def parse_svim(svim_variants):
	sv_dict_total = defaultdict(lambda: defaultdict(list)) # Store all variant coordinates in dictionary using filters and sv types as keys
	sv_dict_total_vcf_line = defaultdict(lambda: defaultdict(list)) # Store all vcf lines in dictionary using filters and sv types as keys
	excluded_due_to_type = []
	excluded_due_to_size = []

	for line in svim_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			line = line.replace(variant_name, "svim_sv")
			variant_qual = int(line_split[5]) # svim qual score
			filter_line = line_split[6] # Describes if variant call passed. Most calls are designated as PASSED, so I probably can ignore the rest, which are probably low quality anyways.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			end_coord = None
			variant_type = None
			variant_size = None

			# Get the end coordinate, variant type, and size
			for x in info_line_split:
				if "END=" in x:
					end_coord = x.split("=")[1]
					variant_size = int(end_coord) - int(start_coord)
				elif "SVTYPE=" in x:
					variant_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = x.split("=")[1]

			if variant_type is not None:
				if variant_type not in SVIM_SV_TYPES:
					excluded_due_to_type.append(line)
				else:
					if variant_type != "BND":
						if variant_type == "INS":
							abs_variant_size= abs(int(variant_size_vcf))
						else:
							variant_size_from_coords = int(end_coord) - int(start_coord)
							abs_variant_size= abs(int(variant_size_from_coords))

						if variant_type in SVIM_SV_TYPES and variant_qual >= args.min_qual_svim and abs_variant_size >= args.minsize:
							if variant_type == "INS":
								bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
								sv_dict_total[filter_line][variant_type].append(bed_line)
								sv_dict_total_vcf_line[filter_line][variant_type].append(line)
							else:
								bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
								sv_dict_total[filter_line][variant_type].append(bed_line)
								sv_dict_total_vcf_line[filter_line][variant_type].append(line)
						else:
							if variant_type not in SVIM_SV_TYPES:
								print("Excluding: " + line)
								print("Invalid variant type: " + variant_type)
							else:
								excluded_due_to_size.append(line)
					else:
						if variant_qual >= args.min_qual_svim:
							end_coord = line_split[4]
							bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
							sv_dict_total[filter_line][variant_type].append(bed_line)
							sv_dict_total_vcf_line[filter_line][variant_type].append(line)
			else:
				print("Variant type is None")
				print(line)

	return (sv_dict_total_vcf_line, sv_dict_total,excluded_due_to_type, excluded_due_to_size)

def get_intersecting_calls_coords(caller_shuffled, caller_original):
	intersection_statistics = [["FILTER,SVTYPE,TOTAL_CALLS,INTERSECTING,NON_INTERSECTING,NON_INTERSECTING_PROPORTION"]] # Store comparison statistics: shuffled_filter, sv_type, total calls, intersecting, non-intersecting, non-intersecting proportion
	non_intersecting_calls = ["readorder,filter,sv type,chr,start,end"] # Store variants found only in original or shuffled
	intersection_statistics_relaxed = [["FILTER,SVTYPE,TOTAL_CALLS,INTERSECTING,NON_INTERSECTING,NON_INTERSECTING_PROPORTION"]]
	non_intersecting_calls_relaxed = ["readorder,filter,sv type,chr,start,end"] # Store variants found only in original or shuffled
	original_filters = caller_original.keys()
	shuffled_filters = caller_shuffled.keys()
	common_filters = set(shuffled_filters) & set(original_filters)

	original_filters_unique = set(original_filters) - set(shuffled_filters)
	for unique_filter in original_filters_unique:
		print("Found filter unique to original vcf file: " + unique_filter)
		for sv_type in caller_original[unique_filter]:
			for line in caller_original[unique_filter][sv_type]:
				non_intersecting_line = "original," + unique_filter + "," + sv_type + "," + str(line).replace("\t", ",")
				non_intersecting_calls.append(non_intersecting_line)

		total_sv_calls_unique_original = len(caller_original[unique_filter][sv_type])
		if total_sv_calls_unique_original > 0:
			intersection_statistics_line = [unique_filter,sv_type, 0, total_sv_calls_unique_original,1]
			intersection_statistics.append(intersection_statistics_line)

	shuffled_filters_unique = set(shuffled_filters) - set(original_filters)
	for unique_filter in shuffled_filters_unique:
		print("Found filter unique to shuffled vcf file: " + unique_filter)
		for sv_type in caller_shuffled[unique_filter]:
			for line in caller_shuffled[unique_filter][sv_type]:
				non_intersecting_line = "shuffled," + unique_filter + "," + sv_type + "," + str(line).replace("\t", ",")
				non_intersecting_calls.append(non_intersecting_line)

		total_sv_calls_unique_shuffled = len(caller_shuffled[unique_filter][sv_type])
		if total_sv_calls_unique_shuffled > 0:
			intersection_statistics_line = [unique_filter,sv_type, 0, total_sv_calls_unique_shuffled,1]
			intersection_statistics.append(intersection_statistics_line)

	for sv_filter in common_filters:
		original_sv_types = caller_original[sv_filter].keys()
		shuffled_sv_types = caller_shuffled[sv_filter].keys()
		common_sv_types = set(original_sv_types) & set(shuffled_sv_types)
		original_sv_types_unique = set(original_sv_types) - set(shuffled_sv_types)
		shuffled_sv_types_unique = set(shuffled_sv_types) - set(original_sv_types)

		for unique_sv in original_sv_types_unique:
			print("Found sv type unique to original vcf file for the " + sv_filter + " filter: " + unique_sv)
			for line in caller_original[sv_filter][unique_sv]:
				non_intersecting_line = "original," + sv_filter + "," + unique_sv + "," + str(line).replace("\t", ",")
				non_intersecting_calls.append(non_intersecting_line)

			total_sv_calls_unique_original = len(caller_original[sv_filter][unique_sv])
			if total_sv_calls_unique_original > 0:
				intersection_statistics_line = [sv_filter,unique_sv, total_sv_calls_unique_original, 0, total_sv_calls_unique_original,1]
				intersection_statistics.append(intersection_statistics_line)

		for unique_sv in shuffled_sv_types_unique:
			print("Found sv type unique to shuffled vcf file for the " + sv_filter + " filter: " + unique_sv)
			for line in caller_shuffled[sv_filter][unique_sv]:
				non_intersecting_line = "shuffled," + sv_filter + "," + unique_sv + "," + str(line).replace("\t", ",")
				non_intersecting_calls.append(non_intersecting_line)

			total_sv_calls_unique_shuffled = len(caller_shuffled[sv_filter][unique_sv])
			if total_sv_calls_unique_shuffled > 0:
				intersection_statistics_line = [sv_filter,unique_sv, total_sv_calls_unique_shuffled, 0, total_sv_calls_unique_shuffled,1]
				intersection_statistics.append(intersection_statistics_line)

		for sv_type in common_sv_types:
				if sv_type != "BND" and sv_type != "INS":
					sv_calls_original = caller_original[sv_filter][sv_type]
					bedlines_original = '\n'.join(str(e) for e in sv_calls_original)
					bedfile_original = BedTool(bedlines_original, from_string=True)

					sv_calls_shuffled = caller_shuffled[sv_filter][sv_type]
					bedlines_shuffled = '\n'.join(str(e) for e in sv_calls_shuffled)
					bedfile_shuffled = BedTool(bedlines_shuffled, from_string=True)
					bedfile_intersection = bedfile_original.intersect(bedfile_shuffled, wa=True, f= 1, r=True)
					bedfile_intersection_no_duplicates = set(bedfile_intersection) # Remove duplicate entries that resulted from multiple overlapping ranges
					bedfile_intersection_relaxed = bedfile_original.intersect(bedfile_shuffled, wa=True, f= 0.5, r=True)
					bedfile_intersection_relaxed_no_duplicates = set(bedfile_intersection_relaxed) # Remove duplicate entries that resulted from multiple overlapping ranges

					original_bedfile_only = bedfile_original.intersect(bedfile_shuffled, wa=True, f= 1, r=True, v = True)
					original_bedfile_only_relaxed = set(original_bedfile_only) - bedfile_intersection_relaxed_no_duplicates

					shuffled_bedfile_only = bedfile_shuffled.intersect(bedfile_original, wa=True, f= 1, r=True, v = True)
					shuffled_bedfile_only_relaxed = set(shuffled_bedfile_only) - bedfile_intersection_relaxed_no_duplicates

					for x in original_bedfile_only:
						non_intersecting_line = "original," + sv_filter + "," + sv_type + "," + str(x).replace("\t", ",")
						non_intersecting_calls.append(non_intersecting_line)

					for x in shuffled_bedfile_only:
						non_intersecting_line = "shuffled," + sv_filter + "," + sv_type + "," + str(x).replace("\t", ",")
						non_intersecting_calls.append(non_intersecting_line)

					for x in original_bedfile_only_relaxed:
						non_intersecting_line = "original," + sv_filter + "," + sv_type + "," + str(x).replace("\t", ",")
						non_intersecting_calls_relaxed.append(non_intersecting_line)

					for x in shuffled_bedfile_only_relaxed:
						non_intersecting_line = "shuffled," + sv_filter + "," + sv_type + "," + str(x).replace("\t", ",")
						non_intersecting_calls_relaxed.append(non_intersecting_line)

					total_sv_calls = len(bedfile_intersection) + len(original_bedfile_only) + len(shuffled_bedfile_only)
					total_sv_calls_relaxed = len(bedfile_intersection_relaxed_no_duplicates) + len(original_bedfile_only_relaxed) + len(shuffled_bedfile_only_relaxed)
					total_non_intersecting = len(original_bedfile_only) + len(shuffled_bedfile_only)
					total_non_intersecting_relaxed = len(original_bedfile_only_relaxed) + len(shuffled_bedfile_only_relaxed)

					if total_sv_calls > 0:
						non_intersecting_proportion = total_non_intersecting / total_sv_calls
						intersection_statistics_line = [sv_filter,sv_type, total_sv_calls, len(bedfile_intersection), total_non_intersecting,non_intersecting_proportion]
						intersection_statistics.append(intersection_statistics_line)
					else:
						non_intersecting_proportion = 0

					if total_sv_calls_relaxed > 0:
						non_intersecting_proportion_relaxed = total_non_intersecting_relaxed / total_sv_calls_relaxed
						intersection_statistics_line = [sv_filter,sv_type, total_sv_calls_relaxed, len(bedfile_intersection_relaxed_no_duplicates), total_non_intersecting_relaxed,non_intersecting_proportion_relaxed]
						intersection_statistics_relaxed.append(intersection_statistics_line)
					else:
						non_intersecting_proportion_relaxed = 0

				else: # This part only applies to exact matches as the relaxed reciprocal overlap analysis does't make sense for single coordinate variants
					original_only = []
					shuffled_only = []
					sv_calls_original = caller_original[sv_filter][sv_type]
					sv_calls_shuffled = caller_shuffled[sv_filter][sv_type]
					sv_calls_original_set = set(sv_calls_original)
					sv_calls_shuffled_set = set(sv_calls_shuffled)
					sv_intersection = sv_calls_original_set.intersection(sv_calls_shuffled_set)
					sv_calls_original_only = sv_calls_original_set.difference(sv_calls_shuffled_set)
					sv_calls_shuffled_only = sv_calls_shuffled_set.difference(sv_calls_original_set)
					total_non_intersecting = len(sv_calls_original_only) + len(sv_calls_shuffled_only)
					total_sv_calls = total_non_intersecting + len(sv_intersection)

					if total_sv_calls > 0:
						non_intersecting_proportion = total_non_intersecting / total_sv_calls
					else:
						non_intersecting_proportion = 0

					intersection_statistics_line = [sv_filter,sv_type, total_sv_calls, len(sv_intersection), total_non_intersecting,non_intersecting_proportion]
					intersection_statistics.append(intersection_statistics_line)

					for line in sv_calls_original_only:
						non_intersecting_line = sv_filter + "," + sv_type + "," + line.replace("\t",",")
						non_intersecting_line = "original," + sv_filter + "," + sv_type + "," + line.replace("\t", ",")
						non_intersecting_calls.append(non_intersecting_line)

					for line in sv_calls_shuffled_only:
						non_intersecting_line = sv_filter + "," + sv_type + "," + line.replace("\t",",")
						non_intersecting_line = "shuffled," + sv_filter + "," + sv_type + "," + line.replace("\t", ",")
						non_intersecting_calls.append(non_intersecting_line)

	total_calls = sum([int(row[2]) for row in intersection_statistics[1:]])
	total_same_breakpoints = sum([int(row[3]) for row in intersection_statistics[1:]])
	total_non_intersecting_calls = sum([int(row[4]) for row in intersection_statistics[1:]])

	if total_calls != 0:
		non_intersecting_proportion = total_non_intersecting_calls / (total_calls)
	else:
		non_intersecting_proportion = 0

	intersection_statistics.append(["Total","ALL_SVS", total_calls, total_same_breakpoints, total_non_intersecting_calls, non_intersecting_proportion])

	total_calls_relaxed = sum([int(row[2]) for row in intersection_statistics_relaxed[1:]])
	total_same_breakpoints_relaxed = sum([int(row[3]) for row in intersection_statistics_relaxed[1:]])
	total_non_intersecting_calls_relaxed = sum([int(row[4]) for row in intersection_statistics_relaxed[1:]])

	if total_calls_relaxed != 0:
		non_intersecting_proportion_relaxed = total_non_intersecting_calls_relaxed / (total_calls_relaxed)
	else:
		non_intersecting_proportion_relaxed = 0

	intersection_statistics_relaxed.append(["Total","ALL_SVS", total_calls_relaxed, total_same_breakpoints_relaxed, total_non_intersecting_calls_relaxed, non_intersecting_proportion_relaxed])

	return intersection_statistics, non_intersecting_calls, intersection_statistics_relaxed, non_intersecting_calls_relaxed

def get_intersecting_calls_vcf_lines(caller_shuffled, caller_original):
	intersection_statistics = [["FILTER,SVTYPE,TOTAL_CALLS,INTERSECTING,NON_INTERSECTING,NON_INTERSECTING_PROPORTION"]] # Store comparison statistics: shuffled_filter, sv_type, total calls, intersecting, non-intersecting, non-intersecting proportion
	non_intersecting_calls = ["readorder,filter,sv type,chr,start,end"] # Store variants found only in original or shuffled
	original_filters = caller_original.keys()
	shuffled_filters = caller_shuffled.keys()
	common_filters = set(shuffled_filters) & set(original_filters)

	original_filters_unique = set(original_filters) - set(shuffled_filters)
	for unique_filter in original_filters_unique:
		print("Found filter unique to original vcf file: " + unique_filter)
		for sv_type in caller_original[unique_filter]:
			for line in caller_original[unique_filter][sv_type]:
				non_intersecting_line = "original," + unique_filter + "," + sv_type + "," + line
				non_intersecting_calls.append(non_intersecting_line)

		total_sv_calls_unique_original = len(caller_original[unique_filter][sv_type])
		if total_sv_calls_unique_original > 0:
			intersection_statistics_line = [sv_filter,sv_type, 0, total_sv_calls_unique_original,1]
			intersection_statistics.append(intersection_statistics_line)

	shuffled_filters_unique = set(shuffled_filters) - set(original_filters)
	for unique_filter in shuffled_filters_unique:
		print("Found filter unique to shuffled vcf file: " + unique_filter)
		for sv_type in caller_shuffled[unique_filter]:
			for line in caller_shuffled[unique_filter][sv_type]:
				non_intersecting_line = "shuffled," + unique_filter + "," + sv_type + "," + line
				non_intersecting_calls.append(non_intersecting_line)

		total_sv_calls_unique_shuffled = len(caller_shuffled[unique_filter][sv_type])
		if total_sv_calls_unique_shuffled > 0:
			intersection_statistics_line = [unique_filter,sv_type, 0, total_sv_calls_unique_shuffled,1]
			intersection_statistics.append(intersection_statistics_line)


	for sv_filter in common_filters:
		original_sv_types = caller_original[sv_filter].keys()
		shuffled_sv_types = caller_shuffled[sv_filter].keys()
		common_sv_types = set(original_sv_types) & set(shuffled_sv_types)
		original_sv_types_unique = set(original_sv_types) - set(shuffled_sv_types)
		shuffled_sv_types_unique = set(shuffled_sv_types) - set(original_sv_types)

		for unique_sv in original_sv_types_unique:
			print("Found sv type unique to original vcf file for the " + sv_filter + " filter: " + unique_sv)
			for line in caller_original[sv_filter][unique_sv]:
				non_intersecting_line = "original," + sv_filter + "," + unique_sv + "," + line
				non_intersecting_calls.append(non_intersecting_line)

			total_sv_calls_unique_original = len(caller_original[sv_filter][unique_sv])
			if total_sv_calls_unique_original > 0:
				intersection_statistics_line = [sv_filter,unique_sv, total_sv_calls_unique_original, 0, total_sv_calls_unique_original,1]
				intersection_statistics.append(intersection_statistics_line)

		for unique_sv in shuffled_sv_types_unique:
			print("Found sv type unique to shuffled vcf file for the " + sv_filter + " filter: " + unique_sv)
			for line in caller_shuffled[sv_filter][unique_sv]:
				non_intersecting_line = "shuffled," + sv_filter + "," + unique_sv + "," + line
				non_intersecting_calls.append(non_intersecting_line)

			total_sv_calls_unique_shuffled = len(caller_shuffled[sv_filter][unique_sv])
			if total_sv_calls_unique_shuffled > 0:
				intersection_statistics_line = [sv_filter,unique_sv, total_sv_calls_unique_shuffled, 0, total_sv_calls_unique_shuffled,1]
				intersection_statistics.append(intersection_statistics_line)

		for sv_type in common_sv_types:
				if sv_type != "BND" and sv_type != "INS":
					sv_calls_original = caller_original[sv_filter][sv_type]
					sv_calls_shuffled = caller_shuffled[sv_filter][sv_type]
					vcf_file_intersection = [line for line in sv_calls_original if line in sv_calls_shuffled]
					vcf_file_intersection_no_duplicates = set(vcf_file_intersection) # Remove duplicate entries that resulted from multiple overlapping ranges

					if len(vcf_file_intersection) != len(vcf_file_intersection_no_duplicates):
						print("Warning intersection between original and shuffled vcf files has duplicates for the filter " + sv_filter + " and sv type " + sv_type)
						print("input vcfs = " + args.original_vcf + " " + args.shuffled_vcf)
						duplicates = [item for item, count in collections.Counter(vcf_file_intersection).items() if count > 1]
						for x in duplicates:
							print(x)

					sv_calls_shuffled_set = set(sv_calls_shuffled)
					original_vcf_file_only = [x for x in sv_calls_original if x not in sv_calls_shuffled_set]
					sv_calls_original_set = set(sv_calls_original)
					shuffled_vcf_file_only = [x for x in sv_calls_shuffled if x not in sv_calls_original_set]

					for x in original_vcf_file_only:
						non_intersecting_line = "original," + sv_filter + "," + sv_type + "," + x
						non_intersecting_calls.append(non_intersecting_line)

					for x in shuffled_vcf_file_only:
						non_intersecting_line = "shuffled," + sv_filter + "," + sv_type + "," + x
						non_intersecting_calls.append(non_intersecting_line)

					total_sv_calls = len(vcf_file_intersection) + len(original_vcf_file_only) + len(shuffled_vcf_file_only)
					total_non_intersecting = len(original_vcf_file_only) + len(shuffled_vcf_file_only)

					if total_sv_calls > 0:
						non_intersecting_proportion = total_non_intersecting / total_sv_calls
						intersection_statistics_line = [sv_filter,sv_type, total_sv_calls, len(vcf_file_intersection), total_non_intersecting,non_intersecting_proportion]
						intersection_statistics.append(intersection_statistics_line)
					else:
						non_intersecting_proportion = 0

				else:
					original_only = []
					shuffled_only = []
					sv_calls_original = caller_original[sv_filter][sv_type]
					sv_calls_shuffled = caller_shuffled[sv_filter][sv_type]
					sv_calls_original_set = set(sv_calls_original)
					sv_calls_shuffled_set = set(sv_calls_shuffled)
					sv_intersection = sv_calls_original_set.intersection(sv_calls_shuffled_set)
					sv_calls_original_only = sv_calls_original_set.difference(sv_calls_shuffled_set)
					sv_calls_shuffled_only = sv_calls_shuffled_set.difference(sv_calls_original_set)
					total_non_intersecting = len(sv_calls_original_only) + len(sv_calls_shuffled_only)
					total_sv_calls = total_non_intersecting + len(sv_intersection)

					if total_sv_calls > 0:
						non_intersecting_proportion = total_non_intersecting / total_sv_calls
					else:
						non_intersecting_proportion = 0

					intersection_statistics_line = [sv_filter,sv_type, total_sv_calls, len(sv_intersection), total_non_intersecting,non_intersecting_proportion]
					intersection_statistics.append(intersection_statistics_line)

					for line in sv_calls_original_only:
						non_intersecting_line = "original," + sv_filter + "," + sv_type + "," + line
						non_intersecting_calls.append(non_intersecting_line)

					for line in sv_calls_shuffled_only:
						non_intersecting_line = "shuffled," + sv_filter + "," + sv_type + "," + line
						non_intersecting_calls.append(non_intersecting_line)

	total_calls = sum([int(row[2]) for row in intersection_statistics[1:]])
	total_same_breakpoints = sum([int(row[3]) for row in intersection_statistics[1:]])
	total_non_intersecting_calls = sum([int(row[4]) for row in intersection_statistics[1:]])

	if total_calls != 0:
		non_intersecting_proportion = total_non_intersecting_calls / (total_calls)
	else:
		non_intersecting_proportion = 0

	intersection_statistics.append(["Total","ALL_SVS", total_calls, total_same_breakpoints, total_non_intersecting_calls, non_intersecting_proportion])

	return intersection_statistics, non_intersecting_calls

# Write summary statistics for overlap analysis
def write_sv_overlap_statistics(intersection_statistics, file_suffix):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/qual_" + str(args.min_qual_svim) + "/summary/"
	else:
		outdir = vcf_dir + "/summary/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/summary_overlap" + file_suffix + ".csv"
	intersection_statistics_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
	header = intersection_statistics[0]
	header = (','.join(str(e) for e in header) + "\n")
	footer = intersection_statistics[-1]
	footer = (','.join(str(e) for e in footer) + "\n")

	for line in intersection_statistics[1:-1]:
		intersection_statistics_string.append(','.join(str(e) for e in line))
	intersection_statistics_sorted = sorted(intersection_statistics_string, key=str.casefold,reverse=True)

	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.write(header)
		f.writelines("%s\n" % l for l in intersection_statistics_sorted)
		f.write(footer)

# Write summary statistics for overlap analysis
def write_breakpoint_statistics(breakpoint_statistics, file_suffix):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/qual_" + str(args.min_qual_svim) + "/summary/"
	else:
		outdir = vcf_dir + "/summary/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/summary_breakpoints" + file_suffix + ".csv"
	breakpoint_statistics_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
	header = breakpoint_statistics[0]
	header = (','.join(str(e) for e in header) + "\n")
	footer = breakpoint_statistics[-1]
	footer = (','.join(str(e) for e in footer) + "\n")

	for line in breakpoint_statistics[1:-1]:
		breakpoint_statistics_string.append(','.join(str(e) for e in line))
	breakpoint_statistics_sorted = sorted(breakpoint_statistics_string, key=str.casefold,reverse=True)

	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.write(header)
		f.writelines("%s\n" % l for l in breakpoint_statistics_sorted)
		f.write(footer)

# Write summary statistics for exact match analysis
def write_non_intersecting_statistics(non_intersecting_statistics, file_suffix):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/qual_" + str(args.min_qual_svim) + "/summary/"
	else:
		outdir = vcf_dir + "/summary/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/overlap_comparison" + file_suffix + ".csv"
	non_intersecting_statistics_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
	header = non_intersecting_statistics[0]
	header = (','.join(str(e) for e in header) + "\n")
	footer = non_intersecting_statistics[-1]
	footer = (','.join(str(e) for e in footer) + "\n")

	for line in non_intersecting_statistics[1:-1]:
		non_intersecting_statistics_string.append(','.join(str(e) for e in line))
	non_intersecting_statistics_sorted = sorted(non_intersecting_statistics_string, key=str.casefold,reverse=True)

	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.write(header)
		f.writelines("%s\n" % l for l in non_intersecting_statistics_sorted)
		f.write(footer)

# Write unique variant calls to disk
def write_unique_calls(unique_calls_original, unique_calls_shuffled, file_suffix):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/qual_" + str(args.min_qual_svim) + "/unique/"
	else:
		outdir = vcf_dir + "/unique/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/original_unique" + file_suffix + ".csv"
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in unique_calls_original)

	output_file = outdir + "/shuffled_unique" + file_suffix + ".csv"
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in unique_calls_shuffled)


# Write unique breakpoints to disk
def write_unique_breakpoints(different_breakpoints, file_suffix):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/qual_" + str(args.min_qual_svim) + "/unique/"
	else:
		outdir = vcf_dir + "/unique/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/different_breakpoints" + file_suffix + ".csv"
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in different_breakpoints)

# Write unique calls to disk
def write_non_intersecting(non_intersecting, file_suffix):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/qual_" + str(args.min_qual_svim) + "/unique/"
	else:
		outdir = vcf_dir + "/unique/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/non_intersecting" + file_suffix + ".csv"
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in non_intersecting)

# Write svs that were excluded to disk
def write_excluded(excluded_lines, filename):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/qual_" + str(args.min_qual_svim) + "/excluded/"
	else:
		outdir = vcf_dir + "/excluded/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/" + filename
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in excluded_lines)

# Function to return summarize the statistics for unique calls
# Unique calls correspond to either:
# sv type not found in calls from the original or shuffled fastqs for a filter found in both of them
# filter not in calls from the original or shuffled fastqs
def get_unique_filter_sv_type_statistics(variant_dict, sv_filter, sv_type):
	unique_calls = [] # Store unique call lines to write to file
	sv_calls_unique = variant_dict[sv_filter][sv_type]

	for line in sv_calls_unique:
		unique_line = sv_filter + "," + sv_type + "," + line.replace("\t",",")
		unique_calls.append(unique_line)

	sv_calls_unique = variant_dict[sv_filter][sv_type]
	if sv_type != "BND":

		bedlines_unique = '\n'.join(str(e) for e in sv_calls_unique)
		bedfile_unique = BedTool(bedlines_unique, from_string=True)
		bedfile_unique_length = len(bedfile_unique)

	else:
		bnd_calls = set(sv_calls_unique)
		bedfile_unique_length = len(bnd_calls)

	return(bedfile_unique_length, unique_calls)

def get_intersection_statistics(caller_original, caller_shuffled, sv_filter, sv_type):
	sv_calls_original = caller_original[sv_filter][sv_type]
	bedlines_original = '\n'.join(str(e) for e in sv_calls_original)
	bedfile_original = BedTool(bedlines_original, from_string=True)
	sv_calls_shuffled = caller_shuffled[sv_filter][sv_type]
	bedlines_shuffled = '\n'.join(str(e) for e in sv_calls_shuffled)
	bedfile_shuffled = BedTool(bedlines_shuffled, from_string=True)

	original_only = []
	shuffled_only = []

	if sv_type != "BND":
		if args.reciprocal is not None:
			bedfile_intersection = bedfile_original.intersect(bedfile_shuffled,wa =True, wb=True, f= args.reciprocal, r=True)
			for x in bedfile_intersection:
				print(x)
				input(".")
			original_only_bedfile = bedfile_original.intersect(bedfile_shuffled, v=True, f= args.reciprocal, r=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
		else:
			bedfile_intersection = bedfile_original.intersect(bedfile_shuffled)
			original_only_bedfile = bedfile_original.intersect(bedfile_shuffled, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.

		for line in original_only_bedfile:
			original_only.append(sv_filter + "," + sv_type + "," + str(line).replace("\t",","))

		if args.reciprocal is not None:
			shuffled_only_bedfile = bedfile_shuffled.intersect(bedfile_original, v=True, f= args.reciprocal, r=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
		else:
			shuffled_only_bedfile = bedfile_shuffled.intersect(bedfile_original, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.

		for line in shuffled_only_bedfile:
			shuffled_only.append(sv_filter + "," + sv_type + "," + str(line).replace("\t",","))
		unique_call_count = len(original_only_bedfile) + len(shuffled_only_bedfile)
		total_call_count = unique_call_count + len(bedfile_intersection)
		unique_call_proportion = unique_call_count / total_call_count
		if unique_call_proportion > 1:
			print("Exceeded 1 for unique call proportion")
		bedfile_intersection_statistics = [sv_filter, sv_type, len(bedfile_intersection), len(original_only_bedfile), len(shuffled_only_bedfile), unique_call_count,unique_call_proportion]

	else:
		bnd_calls_original_set = set(sv_calls_original)
		bnd_calls_shuffled_set = set(sv_calls_shuffled)
		bnd_intersection = bnd_calls_original_set.intersection(bnd_calls_shuffled_set)
		bnd_calls_original_only = bnd_calls_original_set.difference(bnd_calls_shuffled_set)
		bnd_calls_shuffled_only = bnd_calls_shuffled_set.difference(bnd_calls_original_set)
		unique_call_count = len(bnd_calls_original_only) + len(bnd_calls_shuffled_only)
		total_call_count = unique_call_count + len(bnd_intersection)
		unique_call_proportion = unique_call_count / total_call_count
		if unique_call_proportion > 1:
			print("Exceeded 1 for unique call proportion")
		bedfile_intersection_statistics = [sv_filter, sv_type, len(bnd_intersection), len(bnd_calls_original_only ), len(bnd_calls_shuffled_only), unique_call_count, unique_call_proportion]

		for line in bnd_calls_original_only:
			unique_line = sv_filter + "," + sv_type + "," + line.replace("\t",",")
			original_only.append(unique_line)

		for line in bnd_calls_shuffled_only:
			unique_line = sv_filter + "," + sv_type + "," + line.replace("\t",",")
			shuffled_only.append(unique_line)

	return(bedfile_intersection_statistics, original_only, shuffled_only)

def merge_defaultdicts(single_coord_dict,two_coord_dict):
	for vcf_filter in two_coord_dict.keys():
		if (vcf_filter in single_coord_dict):
			for variant_type in two_coord_dict[vcf_filter].keys():
				if (variant_type in single_coord_dict[vcf_filter]):
					print("Same variant type found in both dicts: " + variant_type + " --> dicts should contain different variant types")
				else:
					for bedline in two_coord_dict[vcf_filter][variant_type]:
						single_coord_dict[vcf_filter][variant_type].append(bedline)
		else:
			for variant_type in two_coord_dict[vcf_filter].keys():
				for bedline in two_coord_dict[vcf_filter][variant_type]:
					single_coord_dict[vcf_filter][variant_type].append(bedline)
	return single_coord_dict

with open(args.original_vcf) as f:
	original_variants = f.readlines()
	if args.sv_caller.lower() == "pbsv":
		caller_original = parse_pbsv(original_variants)
	elif args.sv_caller.lower() == "sniffles":
		caller_original = parse_sniffles(original_variants)
	elif args.sv_caller.lower() == "svim":
		caller_original = parse_svim(original_variants)

	excluded_due_to_type = caller_original[2]
	write_excluded(excluded_due_to_type, "original_excluded_due_type.vcf")
	excluded_due_to_size = caller_original[3]
	write_excluded(excluded_due_to_size, "original_excluded_due_size.vcf")

	caller_original_total_svs_vcf_lines = caller_original[0]
	caller_original_total_svs_coords = caller_original[1]

with open(args.shuffled_vcf) as f:
	shuffled_variants = f.readlines()
	if args.sv_caller.lower() == "pbsv":
		caller_shuffled = parse_pbsv(shuffled_variants)
	elif args.sv_caller.lower() == "sniffles":
		caller_shuffled = parse_sniffles(shuffled_variants)
	elif args.sv_caller.lower() == "svim":
		caller_shuffled = parse_svim(shuffled_variants)

	excluded_due_to_type = caller_shuffled[2]
	write_excluded(excluded_due_to_type, "shuffled_excluded_due_type.vcf")
	excluded_due_to_size = caller_shuffled[3]
	write_excluded(excluded_due_to_size, "shuffled_excluded_due_size.vcf")

	caller_shuffled_total_svs_vcf_lines = caller_shuffled[0]
	caller_shuffled_total_svs_coords = caller_shuffled[1]

	intersection_comparison_exact_match = get_intersecting_calls_vcf_lines(caller_shuffled_total_svs_vcf_lines, caller_original_total_svs_vcf_lines)
	intersection_comparison_coords = get_intersecting_calls_coords(caller_shuffled_total_svs_coords, caller_original_total_svs_coords)

	write_non_intersecting_statistics(intersection_comparison_exact_match[0], "_all_svs")
	write_non_intersecting(intersection_comparison_exact_match[1], "_all_svs")
	write_non_intersecting_statistics(intersection_comparison_coords[0], "_coords_all_svs")
	write_non_intersecting(intersection_comparison_coords[1], "_coords_all_svs")
	write_non_intersecting_statistics(intersection_comparison_coords[2], "_relaxed_all_svs")
	write_non_intersecting(intersection_comparison_coords[3], "_relaxed_all_svs")
