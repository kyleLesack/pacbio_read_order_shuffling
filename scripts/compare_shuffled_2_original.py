import argparse
import os
from pybedtools import BedTool

parser = argparse.ArgumentParser()
parser.add_argument("shuffled_vcf", help="vcf file created from shuffled FASTQ")
parser.add_argument("original_vcf", help="vcf file created from original FASTQ") # The original FASTQ was used to create shuffled FASTQ files
parser.add_argument("sv_caller", help="variant caller used to create vcf file. Use pbsv, sniffles, or svim")
parser.add_argument("--minsize", default=0, type=int, help="minimum variant size in bp")
parser.add_argument("--maxsize", default=500000, type=int, help="maximum variant size in bp")
parser.add_argument("--min_qual_svim", default=0, type=int, help="minimum qual value for svim")
args = parser.parse_args()

# Types of variants to look for
PBSV_SV_TYPES = ["INV","DEL","DUP", "INS", "BND"]
sniffles_SV_TYPES = ["INV", "DEL", "DUP", "INVDUP", "INV/INVDUP", "DEL/INV", "TRA", "INS", "BND", "DUP/INS"] # Note: DUP/INS causes problems, as some have END coords that start before the start. DUP/INS is rare, so I can ignore them altogether.
SVIM_SV_TYPES = ["INV", "DEL","DUP:INT", "DUP:TANDEM","INS","DUP","BND"]

from collections import defaultdict

# Extract pbsv variants
def parse_pbsv(pbsv_variants):
	callerDict = defaultdict(lambda: defaultdict(list)) # Store variants in dict

	for line in pbsv_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			#variant_qual = int(line_split[5]) # pbsv qual score
			filter_line = line_split[6] # Describes if variant call passed.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			if info_line_split[0] == "IMPRECISE":
				filter_line = filter_line + "-IMPRECISE"
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
				if variant_type != "BND":
					if variant_type == "INS":
						abs_variant_size= abs(int(variant_size_vcf))
					else:
						variant_size_from_coords = int(end_coord) - int(start_coord)
						abs_variant_size= abs(int(variant_size_from_coords))

					if variant_type in PBSV_SV_TYPES and abs_variant_size >= args.minsize:
						#bed_line = chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_name + "\tSupport:" + str(variant_qual)
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][variant_type].append(bed_line)
					else:
						print("Excluding: " + line)
						if variant_type not in pbsv_SV_TYPES:
							print("Invalid variant type: " + variant_type)
						#if variant_qual < args.min_qual_pbsv:
						#	print("pbsv Qual too low: " + variant_qual)
						if abs_variant_size_from_coords < args.minsize:
							print("Variant too small: " + str(abs_variant_size_from_coords))
				else:
					#if variant_qual >= args.min_qual_pbsv:
					end_coord = line_split[4]
					#bed_line = chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_name + "\tSupport:" + str(variant_qual)
					bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
					callerDict[filter_line][variant_type].append(bed_line)
			else:
				print("Variant type is none")
				print(line)

	return (callerDict)

# Extract variant calls from sniffles
def parse_sniffles(sniffles_variants):
	callerDict = defaultdict(lambda: defaultdict(list)) # Store calls in dictionary
	excluded_sv_types = set([])

	for line in sniffles_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			filter_line = line_split[6] # Describes if variant call passed.
			info_line = line_split[7]
			info_line_split = info_line.split(";")
			variant_precision = info_line_split[0] # Describes if breakpoints are precise
			filter_line = filter_line + "-" + variant_precision # Combine filter with precision

			end_coord = None
			sv_type = None
			variant_size_vcf = None

			# Get the end coordinate, variant type, and size
			for x in info_line_split:
				if "END=" in x:
					end_coord = x.split("=")[1]
				elif "SVTYPE=" in x:
					sv_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = x.split("=")[1]

			variant_renamed = "sniffles_" + sv_type + "_" + variant_name
			support_line = line_split[9] # Line that provides read support for the reference allele and sv.
			sv_read_support = support_line.split(":")[2]

			if sv_type=="BND":
				end_coord = line_split[4]
				if end_coord is not None:
					if sv_type in sniffles_SV_TYPES:
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][sv_type].append(bed_line)
					else:
						print("Excluding variant type: " + sv_type + " in vcf file: " + args.input_vcf)
				else:
					print("Excluding: " + line)
					if end_coord is None:
						print("End coordinate is None")

			elif sv_type is not None and end_coord is not None and variant_size_vcf is not None:
				if sv_type in sniffles_SV_TYPES:
					abs_variant_size= abs(int(variant_size_vcf))

					if abs_variant_size >= args.minsize: # Select variants of a minimum size
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][sv_type].append(bed_line)
					else:
						print("Excluding due to size: " + line)
				else:
					print("Excluding variant type: " + sv_type + " in vcf file: " + args.input_vcf)
			else:
				print("Excluding: " + line)
				if sv_type is None:
					print("Variant type is None")
				if end_coord is None:
					print("End coordinate is None")
				if variant_size_vcf is None:
					print("Variant size is None")

	return (callerDict)

# Extract svim variants
def parse_svim(svim_variants):
	callerDict = defaultdict(lambda: defaultdict(list)) # Store variants in dict
	for line in svim_variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			variant_qual = int(line_split[5]) # svim qual score
			filter_line = line_split[6] # Describes if variant call passed. Most calls are designated as PASSED, so I probably can ignore the rest, which are probably low quality anyways.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			end_coord = None
			sv_type = None
			variant_size = None

			# Get the end coordinate, variant type, and size
			for x in info_line_split:
				if "END=" in x:
					end_coord = x.split("=")[1]
					variant_size = int(end_coord) - int(start_coord)
				elif "SVTYPE=" in x:
					sv_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = x.split("=")[1]

			if sv_type is not None:
				if sv_type != "BND":
					if sv_type == "INS":
						abs_variant_size= abs(int(variant_size_vcf))
					else:
						variant_size_from_coords = int(end_coord) - int(start_coord)
						abs_variant_size= abs(int(variant_size_from_coords))

					if sv_type in SVIM_SV_TYPES and variant_qual >= 15 and abs_variant_size >= 0:
						#bed_line = chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_name + "\tSupport:" + str(variant_qual)
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][sv_type].append(bed_line)
					else:
						#print("Excluding: " + line)
						if sv_type not in SVIM_SV_TYPES:
							print("Excluding: " + line)
							print("Invalid variant type: " + sv_type)
						#if variant_qual < args.min_qual_svim:
						#	print("SVIM Qual too low: " + str(variant_qual))
						#if abs_variant_size < args.minsize:
						#	print("Variant too small: " + str(abs_variant_size_from_coords))
				else:
					if variant_qual >= 15:
						end_coord = line_split[4]
						#bed_line = chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_name + "\tSupport:" + str(variant_qual)
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][sv_type].append(bed_line)
			else:
				print("Variant type is None")
				print(line)

	return (callerDict)

def compare_vcf_breakpoints(caller_shuffled, caller_original):
	intersection_statistics = [["FILTER,SVTYPE,INTERSECTION,SAME_BREAKPOINTS,DIFFERENT_BREAKPOINTS,DISCORDANT_PROPORTION"]] # Store comparison statistics: shuffled_filter, sv_type, intersection size, unique to original, unique to shuffled
	different_breakpoints = [] # Store variants called by both the original and shuffled VCF files but with different breakpoints
	original_filters = caller_original.keys()
	shuffled_filters = caller_shuffled.keys()
	#unique_filters = set(shuffled_filters) ^ set(original_filters)
	common_filters = set(shuffled_filters) & set(original_filters)
	for sv_filter in common_filters:
		original_sv_types = caller_original[sv_filter].keys()
		shuffled_sv_types = caller_shuffled[sv_filter].keys()
		common_sv_types = set(original_sv_types) & set(shuffled_sv_types)
		for sv_type in common_sv_types:
				if sv_type != "BND" and sv_type != "INS":
					same_breakpoint_count = 0
					different_breakpoint_count = 0

					sv_calls_original = caller_original[sv_filter][sv_type]
					bedlines_original = '\n'.join(str(e) for e in sv_calls_original)
					bedfile_original = BedTool(bedlines_original, from_string=True)

					sv_calls_shuffled = caller_shuffled[sv_filter][sv_type]
					bedlines_shuffled = '\n'.join(str(e) for e in sv_calls_shuffled)
					bedfile_shuffled = BedTool(bedlines_shuffled, from_string=True)

					bedfile_intersection = bedfile_original.intersect(bedfile_shuffled, wa=True, wb=True)

					for x in bedfile_intersection:
						#chromosome = x[0]
						original_start = x[1]
						original_end = x[2]
						shuffled_start = x[4]
						shuffled_end = x[5]
						if original_start != shuffled_start or original_end != shuffled_end:
							different_breakpoint_count += 1
							different_breakpoints.append(x)
						else:
							same_breakpoint_count +=1
					if len(bedfile_intersection) > 0:
						discordant_proportion = different_breakpoint_count / len(bedfile_intersection)
						intersection_statistics_line = [sv_filter,sv_type, len(bedfile_intersection),same_breakpoint_count,different_breakpoint_count,discordant_proportion]
						intersection_statistics.append(intersection_statistics_line)
					#intersection_statistics.append(','.join(str(e) for e in intersection_statistics_line))
	total_intersection = sum([int(row[2]) for row in intersection_statistics[1:]])
	total_same_breakpoints = sum([int(row[3]) for row in intersection_statistics[1:]])
	total_different_breakpoints = sum([int(row[4]) for row in intersection_statistics[1:]])
	total_discordant = total_different_breakpoints / total_intersection
	intersection_statistics.append(["Total","ALL_SVS", total_intersection, total_same_breakpoints, total_different_breakpoints, total_discordant])

	return intersection_statistics, different_breakpoints

# Write summary statistics to disk
def write_result_statistics(intersection_statistics_total_calls,intersection_statistics_breakpoints):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/QUAL_" + str(args.min_qual_svim) + "/summary/"
	else:
		outdir = vcf_dir + "/summary/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/summary_total_svs.csv"
	intersection_statistics_total_calls_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
	header = intersection_statistics_total_calls[0]
	header = (','.join(str(e) for e in header) + "\n")
	footer = intersection_statistics_total_calls[-1]
	footer = (','.join(str(e) for e in footer) + "\n")

	for line in intersection_statistics_total_calls[1:-1]:
		intersection_statistics_total_calls_string.append(','.join(str(e) for e in line))
	intersection_statistics_sorted = sorted(intersection_statistics_total_calls_string, key=str.casefold,reverse=True)

	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.write(header)
		f.writelines("%s\n" % l for l in intersection_statistics_sorted)
		f.write(footer)

	output_file = outdir + "/summary_breakpoints.csv"
	intersection_statistics_breakpoints_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
	header = intersection_statistics_breakpoints[0]
	header = (','.join(str(e) for e in header) + "\n")
	footer = intersection_statistics_breakpoints[-1]
	footer = (','.join(str(e) for e in footer) + "\n")

	for line in intersection_statistics_breakpoints[1:-1]:
		intersection_statistics_breakpoints_string.append(','.join(str(e) for e in line))
	intersection_statistics_breakpoints_sorted = sorted(intersection_statistics_breakpoints_string, key=str.casefold,reverse=True)

	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.write(header)
		f.writelines("%s\n" % l for l in intersection_statistics_breakpoints_sorted)
		f.write(footer)

# Write variant calls to disk
def write_unique(unique_calls_original, unique_calls_shuffled, different_breakpoints):
	vcf_dir = os.path.dirname(args.shuffled_vcf) # use vcf filename for bedfile
	if args.sv_caller.lower() == "svim":
		outdir = vcf_dir + "/QUAL_" + str(args.min_qual_svim) + "/unique/"
	else:
		outdir = vcf_dir + "/unique/"

	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/original_unique.csv"
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in unique_calls_original)

	output_file = outdir + "/shuffled_unique.csv"
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in unique_calls_shuffled)

	#unique_calls_original_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
	#header = intersection_statistics_total_calls[0]
	#header = (','.join(str(e) for e in header) + "\n")
	#footer = intersection_statistics_total_calls[-1]
	#footer = (','.join(str(e) for e in footer) + "\n")

	#for line in intersection_statistics_total_calls[1:-1]:
	#	intersection_statistics_total_calls_string.append(','.join(str(e) for e in line))
	#intersection_statistics_sorted = sorted(intersection_statistics_total_calls_string, key=str.casefold,reverse=True)


	#with open(output_file, 'w', newline='\n') as f:
	#	f.write(header)
	#	f.writelines("%s\n" % l for l in intersection_statistics_sorted)
	#	f.write(footer)


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
		bedfile_intersection = bedfile_original.intersect(bedfile_shuffled)
		original_only_bedfile = bedfile_original.intersect(bedfile_shuffled, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
		for line in original_only_bedfile:
			original_only.append(sv_filter + "," + sv_type + "," + str(line).replace("\t",","))

		shuffled_only_bedfile = bedfile_shuffled.intersect(bedfile_original, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
		for line in shuffled_only_bedfile:
			shuffled_only.append(sv_filter + "," + sv_type + "," + str(line).replace("\t",","))
		unique_call_count = len(original_only_bedfile) + len(shuffled_only_bedfile)
		total_call_count = unique_call_count + len(bedfile_intersection)
		unique_call_proportion = unique_call_count / total_call_count
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
		bedfile_intersection_statistics = [sv_filter, sv_type, len(bnd_intersection), len(bnd_calls_original_only ), len(bnd_calls_shuffled_only), unique_call_count, unique_call_proportion]

		for line in bnd_calls_original_only:
			unique_line = sv_filter + "," + sv_type + "," + line.replace("\t",",")
			original_only.append(unique_line)

		for line in bnd_calls_shuffled_only:
			unique_line = sv_filter + "," + sv_type + "," + line.replace("\t",",")
			shuffled_only.append(unique_line)

	return(bedfile_intersection_statistics, original_only, shuffled_only)

def compare_vcf_total_predictions_new(caller_shuffled, caller_original):
	#print("Length caller_shuffled (1): " + str(len(caller_shuffled)))
	warning_messages = set()
	unique_calls_original = [] # Store calls unique to the original vcf in list, excludes BND, INS
	unique_calls_shuffled = [] # Store calls unique to the shuffled vcf in list, excludes BND, INS
	intersection_statistics = [["FILTER","SVTYPE","INTERSECTION","ORIGINAL_ONLY","SHUFFLED_ONLY","UNIQUE","UNIQUE_PROPORTION"]] # Store comparison statistics: shuffled_filter, sv_type, intersection size, unique to original, unique to shuffled
	shuffled_filters = caller_shuffled.keys()
	original_filters = caller_original.keys()
	common_filters = set(original_filters) & set(shuffled_filters)
	filter_original_only = set(original_filters) - set(shuffled_filters)
	filter_shuffled_only = set(shuffled_filters) - set(original_filters)
	#print("Length caller_shuffled (2): " + str(len(caller_shuffled)))

	if len(filter_original_only) > 0: # Find variant calls with filters only found using the original fastq file
		for original_filter in filter_original_only:
			warning_messages.add("Filter: " + original_filter + " not found in shuffled vcf file")
			for sv_type in caller_original[original_filter].keys():
				unique_calls_original_stats = get_unique_filter_sv_type_statistics(caller_original, original_filter, sv_type)
				bedfile_intersection_statistics = [original_filter, sv_type, 0, unique_calls_original_stats[0], 0, unique_calls_original_stats[0], 1]
				#print(bedfile_intersection_statistics)
				intersection_statistics.append(bedfile_intersection_statistics)
				unique_calls_original.extend(unique_calls_original_stats[1])

	if len(filter_shuffled_only) > 0: # Find variant calls with filters only found using the shuffled fastq file
		for shuffled_filter in filter_shuffled_only:
			warning_messages.add("Filter: " + shuffled_filter + " not found in original vcf file")
			for sv_type in caller_shuffled[shuffled_filter].keys():
				unique_calls_shuffled_stats = get_unique_filter_sv_type_statistics(caller_shuffled, shuffled_filter, sv_type)
				bedfile_intersection_statistics = [shuffled_filter, sv_type, 0, 0, unique_calls_shuffled_stats[0],unique_calls_shuffled_stats[0], 1]
				#print(bedfile_intersection_statistics)
				intersection_statistics.append(bedfile_intersection_statistics)
				unique_calls_shuffled.extend(unique_calls_shuffled_stats[1])

	for common_filter in common_filters: # Find variant calls for shared filters
		original_sv_types = caller_original[common_filter].keys()
		shuffled_sv_types = caller_shuffled[common_filter].keys()
		sv_original_only = set(original_sv_types) - set(shuffled_sv_types)
		sv_shuffled_only = set(shuffled_sv_types) - set(original_sv_types)
		common_svs = set(original_sv_types) & set(shuffled_sv_types)

		for sv_type in sv_original_only:
			warning_messages.add("SV Type: " + sv_type + " not found in shuffled vcf file")
			unique_calls_original_stats = get_unique_filter_sv_type_statistics(caller_original, common_filter, sv_type)
			bedfile_intersection_statistics = [common_filter, sv_type, 0, unique_calls_original_stats[0], 0,unique_calls_original_stats[0], 1]
			#print(bedfile_intersection_statistics)
			intersection_statistics.append(bedfile_intersection_statistics)

			unique_calls_original.extend(unique_calls_original_stats[1])

		for sv_type in sv_shuffled_only:
			warning_messages.add("SV Type: " + sv_type + " not found in original vcf file")
			unique_calls_shuffled_stats = get_unique_filter_sv_type_statistics(caller_shuffled, common_filter, sv_type)
			bedfile_intersection_statistics = [common_filter, sv_type, 0, 0, unique_calls_shuffled_stats[0],unique_calls_shuffled_stats[0], 1]
			#print(bedfile_intersection_statistics)
			intersection_statistics.append(bedfile_intersection_statistics)
			unique_calls_shuffled.extend(unique_calls_shuffled_stats[1])

		for sv_type in common_svs:
			bedfile_intersection = get_intersection_statistics(caller_original, caller_shuffled, common_filter, sv_type)
			#print("Here common")
			#print(bedfile_intersection[0])
			intersection_statistics.append(bedfile_intersection[0])
			if len(bedfile_intersection[1]) != 0:
				unique_calls_original.extend(bedfile_intersection[1])
			if len(bedfile_intersection[2]) != 0:
				unique_calls_shuffled.extend(bedfile_intersection[2])

	total_intersection = sum([int(row[2]) for row in intersection_statistics[1:]])
	total_original_only = sum([int(row[3]) for row in intersection_statistics[1:]])
	total_shuffled_only = sum([int(row[4]) for row in intersection_statistics[1:]])
	total_unique = sum([int(row[5]) for row in intersection_statistics[1:]])
	total_unique_proportion = total_unique / total_intersection
	intersection_statistics.append(["Total","ALL_SVS", total_intersection, total_original_only, total_shuffled_only, total_unique, total_unique_proportion])

	return(intersection_statistics, unique_calls_original, unique_calls_shuffled)

if args.sv_caller.lower() == "pbsv":
	with open(args.original_vcf) as f:
		original_variants = f.readlines()
		pbsv_original = parse_pbsv(original_variants)

	with open(args.shuffled_vcf) as f:
		shuffled_variants = f.readlines()
		pbsv_shuffled = parse_pbsv(shuffled_variants)

		vcf_comparison_total_svs = compare_vcf_total_predictions_new(pbsv_shuffled, pbsv_original)
		vcf_comparison_breakpoints = compare_vcf_breakpoints(pbsv_shuffled, pbsv_original)

		write_result_statistics(vcf_comparison_total_svs[0], vcf_comparison_breakpoints[0])
		write_unique(vcf_comparison_total_svs[1], vcf_comparison_total_svs[2], vcf_comparison_breakpoints[1])

elif args.sv_caller.lower() == "svim":
	with open(args.original_vcf) as f:
		original_variants = f.readlines()
		svim_original = parse_svim(original_variants)

	with open(args.shuffled_vcf) as f:
		shuffled_variants = f.readlines()
		svim_shuffled = parse_svim(shuffled_variants)

		vcf_comparison_total_svs = compare_vcf_total_predictions_new(svim_shuffled, svim_original)
		vcf_comparison_breakpoints = compare_vcf_breakpoints(svim_shuffled, svim_original)

		write_result_statistics(vcf_comparison_total_svs[0], vcf_comparison_breakpoints[0])
		write_unique(vcf_comparison_total_svs[1], vcf_comparison_total_svs[2], vcf_comparison_breakpoints[1])

elif args.sv_caller.lower() == "sniffles":
	with open(args.original_vcf) as f:
		original_variants = f.readlines()
		sniffles_original = parse_sniffles(original_variants)

	with open(args.shuffled_vcf) as f:
		shuffled_variants = f.readlines()
		sniffles_shuffled = parse_sniffles(shuffled_variants)

		vcf_comparison_total_svs = compare_vcf_total_predictions_new(sniffles_shuffled, sniffles_original)
		vcf_comparison_breakpoints = compare_vcf_breakpoints(sniffles_shuffled, sniffles_original)

		write_result_statistics(vcf_comparison_total_svs[0], vcf_comparison_breakpoints[0])
		write_unique(vcf_comparison_total_svs[1], vcf_comparison_total_svs[2], vcf_comparison_breakpoints[1])
