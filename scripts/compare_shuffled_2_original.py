import argparse
import os
import pybedtools
from pybedtools import BedTool

parser = argparse.ArgumentParser()
parser.add_argument("shuffled_vcf", help="vcf file created from shuffled FASTQ")
parser.add_argument("original_vcf", help="vcf file created from original FASTQ") # The original FASTQ was used to create shuffled FASTQ files
parser.add_argument("sv_caller", help="variant caller used to create vcf file. Use pbsv, sniffles, or svim")
parser.add_argument("--minsize", default=0, type=int, help="minimum variant size in bp")
#parser.add_argument("--maxsize", default=500000, type=int, help="maximum variant size in bp")
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

	for line in variants:
		if line[0] != "#":
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			variant_name = line_split[2]
			#variant_qual = int(line_split[5]) # pbsv qual score
			filter_line = line_split[6] # Describes if variant call passed. Most calls are designated as PASSED, so I probably can ignore the rest, which are probably low quality anyways.
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
	excluded_variant_types = set([])

	for line in variants:
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
			variant_type = None
			variant_size_vcf = None

			# Get the end coordinate, variant type, and size
			for x in info_line_split:
				if "END=" in x:
					end_coord = x.split("=")[1]
				elif "SVTYPE=" in x:
					variant_type = x.split("=")[1]
				elif "SVLEN=" in x:
					variant_size_vcf = x.split("=")[1]

			variant_renamed = "sniffles_" + variant_type + "_" + variant_name
			support_line = line_split[9] # Line that provides read support for the reference allele and sv.
			sv_read_support = support_line.split(":")[2]

			if variant_type=="BND":
				end_coord = line_split[4]
				if end_coord is not None:
					if variant_type in sniffles_SV_TYPES:
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][variant_type].append(bed_line)
					else:
						print("Excluding variant type: " + variant_type + " in vcf file: " + args.input_vcf)
				else:
					print("Excluding: " + line)
					if end_coord is None:
						print("End coordinate is None")

			elif variant_type is not None and end_coord is not None and variant_size_vcf is not None:
				if variant_type in sniffles_SV_TYPES:
					abs_variant_size= abs(int(variant_size_vcf))

					if abs_variant_size >= args.minsize: # Select variants of a minimum size
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][variant_type].append(bed_line)
					else:
						print("Excluding due to size: " + line)
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

	return (callerDict)

# Extract svim variants
def parse_svim(svim_variants):
	callerDict = defaultdict(lambda: defaultdict(list)) # Store variants in dict

	for line in variants:
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

					if variant_type in SVIM_SV_TYPES and variant_qual >= args.min_qual_svim and abs_variant_size >= args.minsize:
						#bed_line = chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_name + "\tSupport:" + str(variant_qual)
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][variant_type].append(bed_line)
					else:
						#print("Excluding: " + line)
						if variant_type not in SVIM_SV_TYPES:
							print("Excluding: " + line)
							print("Invalid variant type: " + variant_type)
						#if variant_qual < args.min_qual_svim:
						#	print("SVIM Qual too low: " + str(variant_qual))
						#if abs_variant_size < args.minsize:
						#	print("Variant too small: " + str(abs_variant_size_from_coords))
				else:
					if variant_qual >= args.min_qual_svim:
						end_coord = line_split[4]
						#bed_line = chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + variant_name + "\tSupport:" + str(variant_qual)
						bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
						callerDict[filter_line][variant_type].append(bed_line)
			else:
				print("Variant type is None")
				print(line)

	return (callerDict)

def compare_vcf_total_predictions(svim_shuffled, svim_original):
	warning_messages = set()
	unique_calls_original = [] # Store calls unique to the original vcf in list, excludes BND, INS
	unique_calls_shuffled = [] # Store calls unique to the shuffled vcf in list, excludes BND, INS
	intersection_statistics = [["FILTER","SVTYPE","INTERSECTION","ORIGINAL_ONLY","SHUFFLED_ONLY"]] # Store comparison statistics: shuffled_filter, variant_type, intersection size, unique to original, unique to shuffled
	shuffled_filters = svim_shuffled.keys()
	original_filters = svim_original.keys()
	unique_filters = set(shuffled_filters) ^ set(original_filters)
	filter_original_only = set(original_filters) - set(shuffled_filters)

	if len(unique_filters) != 0:
		for filter in unique_filters:
			if filter not in shuffled_filters:
				warning_messages.add(("Filter: " + filter + " not found in the shuffled vcf file"))
			if filter not in original_filters:
				warning_messages.add(("Filter: " + filter + " not found in the original vcf file"))

		if len(filter_original_only) > 0:
			for original_filter in svim_original.keys():
				original_variant_types = svim_original[original_filter].keys()
				shuffled_variant_types = svim_shuffled[original_filter].keys()
				unique_variant_types = set(original_variant_types) ^ set(shuffled_variant_types) # Identify any variant types only found in the shuffled or original vcf files

				if len(unique_variant_types) != 0:
					for svtype in unique_variant_types:
						if svtype not in original_variant_types:
							warning_messages.add(("Variant type: " + svtype + " not found in the original vcf file for the filter type: " + original_filter))
						if svtype not in shuffled_variant_types:
							warning_messages.add(("Variant type: " + svtype + " not found in the shuffled vcf file for the filter type: " + original_filter))

				for variant_type in original_variant_types:
					if variant_type != "BND":
						sv_calls_original = svim_original[original_filter][variant_type]
						bedlines_original = '\n'.join(str(e) for e in sv_calls_original)
						bedfile_original = BedTool(bedlines_original, from_string=True)

						if original_filter not in shuffled_filters:
							original_only = bedfile_original
							bedfile_intersection_statistics = [original_filter, variant_type, 0, len(original_only), 0]
							shuffled_only = None
							#print(bedfile_intersection_statistics)
							intersection_statistics.append(bedfile_intersection_statistics)

						elif variant_type not in shuffled_variant_types:
							original_only = bedfile_original
							bedfile_intersection_statistics = [original_filter, variant_type, 0, len(original_only), 0]
							shuffled_only = None
							#print(bedfile_intersection_statistics)
							intersection_statistics.append(bedfile_intersection_statistics)

						else:
							sv_calls_shuffled = svim_shuffled[original_filter][variant_type]
							bedlines_shuffled = '\n'.join(str(e) for e in sv_calls_shuffled)
							bedfile_shuffled = BedTool(bedlines_shuffled, from_string=True)
							bedfile_intersection = bedfile_shuffled.intersect(bedfile_original)
							shuffled_only = bedfile_shuffled.intersect(bedfile_original, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
							original_only = bedfile_original.intersect(bedfile_shuffled, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
							bedfile_intersection_statistics = [original_filter, variant_type, len(bedfile_intersection), len(shuffled_only), len(original_only)]
							#print(bedfile_intersection_statistics)
							intersection_statistics.append(bedfile_intersection_statistics)

						if shuffled_only is not None  and len(shuffled_only) > 0:
							for x in shuffled_only:
								unique_entry = ["shuffled VCF", original_filter, variant_type, x.chrom, x.start, x.stop, len(x)]
								unique_calls_shuffled.append(','.join(str(e) for e in unique_entry))

						if len(original_only) > 0:
							for x in original_only:
								unique_entry = ["original VCF", original_filter, variant_type, x.chrom, x.start, x.stop, len(x)]
								unique_calls_original.append(','.join(str(e) for e in unique_entry))
					else:
						bnd_calls_original = svim_original[original_filter][variant_type]
						bnd_calls_original_set = set(bnd_calls_original)

						if original_filter not in shuffled_filters:
							bnd_calls_original_only = bnd_calls_original
							bedfile_intersection_statistics = [original_filter, variant_type, 0, len(bnd_calls_original_only), 0]
							bnd_calls_shuffled_only = None
						elif variant_type not in shuffled_variant_types:
							bnd_calls_original_only = bnd_calls_original
							bedfile_intersection_statistics = [original_filter, variant_type, 0, len(bnd_calls_original_only), 0]
							bnd_calls_shuffled_only = None

						else:
							bnd_calls_shuffled = svim_shuffled[original_filter][variant_type]
							bnd_calls_shuffled_set = set(bnd_calls_shuffled)
							bnd_intersection = bnd_calls_shuffled_set.intersection(bnd_calls_original_set)
							bnd_calls_shuffled_only = bnd_calls_shuffled_set.difference(bnd_calls_original_set)
							bnd_calls_original_only = bnd_calls_original_set.difference(bnd_calls_shuffled_set)

							bedfile_intersection_statistics = [original_filter, variant_type, len(bnd_intersection), len(bnd_calls_shuffled_only ), len(bnd_calls_original_only)]
							intersection_statistics.append(bedfile_intersection_statistics)

						if len(bnd_calls_shuffled_only) > 0:
							for x in bnd_calls_shuffled_only:
								x_split = x.split("\t")
								unique_entry = ["shuffled VCF", original_filter, variant_type, x_split[0], x_split[1], x_split[2], "NA"]
								unique_calls_shuffled.append(','.join(str(e) for e in unique_entry))

						if bnd_calls_shuffled_only is not None	and	len(bnd_calls_original_only) > 0:
							for x in bnd_calls_original_only:
								x_split = x.split("\t")
								unique_entry = ["original VCF", original_filter, variant_type, x_split[0], x_split[1], x_split[2], "NA"]
								unique_calls_original.append(','.join(str(e) for e in unique_entry))

		for shuffled_filter in svim_shuffled.keys():
			shuffled_variant_types = svim_shuffled[shuffled_filter].keys()
			original_variant_types = svim_original[shuffled_filter].keys()
			unique_variant_types = set(shuffled_variant_types) ^ set(original_variant_types) # Identify any variant types only found in the original or shuffled vcf files

			# NOTE: Add code that still get intersection statistics when the variant types aren't the same
			if len(unique_variant_types) != 0:
				for svtype in unique_variant_types:
					if svtype not in original_variant_types:
						warning_messages.add(("Variant type: " + svtype + " not found in the original vcf file for the filter type: " + shuffled_filter))
					if svtype not in shuffled_variant_types:
						warning_messages.add(("Variant type: " + svtype + " not found in the shuffled vcf file for the filter type: " + shuffled_filter))

			for variant_type in shuffled_variant_types:
				if variant_type != "BND":
					sv_calls_shuffled = svim_shuffled[shuffled_filter][variant_type]
					bedlines_shuffled = '\n'.join(str(e) for e in sv_calls_shuffled)
					bedfile_shuffled = BedTool(bedlines_shuffled, from_string=True)

					if shuffled_filter not in original_filters:
						shuffled_only = bedfile_shuffled
						bedfile_intersection_statistics = [shuffled_filter, variant_type, 0, 0, len(shuffled_only)]
						original_only = None
						#print(bedfile_intersection_statistics)
						intersection_statistics.append(bedfile_intersection_statistics)

					elif variant_type not in original_variant_types:
						shuffled_only = bedfile_shuffled
						bedfile_intersection_statistics = [shuffled_filter, variant_type, 0, 0, len(shuffled_only)]
						original_only = None
						#print(bedfile_intersection_statistics)
						intersection_statistics.append(bedfile_intersection_statistics)

					else:
						sv_calls_original = svim_original[shuffled_filter][variant_type]
						bedlines_original = '\n'.join(str(e) for e in sv_calls_original)
						bedfile_original = BedTool(bedlines_original, from_string=True)
						bedfile_intersection = bedfile_original.intersect(bedfile_shuffled)
						original_only = bedfile_original.intersect(bedfile_shuffled, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
						shuffled_only = bedfile_shuffled.intersect(bedfile_original, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
						bedfile_intersection_statistics = [shuffled_filter, variant_type, len(bedfile_intersection), len(original_only), len(shuffled_only)]
						#print(bedfile_intersection_statistics)
						intersection_statistics.append(bedfile_intersection_statistics)

					if original_only is not None  and len(original_only) > 0:
						for x in original_only:
							unique_entry = ["Original VCF", shuffled_filter, variant_type, x.chrom, x.start, x.stop, len(x)]
							unique_calls_original.append(','.join(str(e) for e in unique_entry))

					if len(shuffled_only) > 0:
						for x in shuffled_only:
							unique_entry = ["Shuffled VCF", shuffled_filter, variant_type, x.chrom, x.start, x.stop, len(x)]
							unique_calls_shuffled.append(','.join(str(e) for e in unique_entry))
				else:
					bnd_calls_shuffled = svim_shuffled[shuffled_filter][variant_type]
					bnd_calls_shuffled_set = set(bnd_calls_shuffled)

					if shuffled_filter not in original_filters:
						bnd_calls_shuffled_only = bnd_calls_shuffled
						bedfile_intersection_statistics = [shuffled_filter, variant_type, 0, 0, len(bnd_calls_shuffled_only)]
						bnd_calls_original_only = None
					else:
						bnd_calls_original = svim_original[shuffled_filter][variant_type]
						bnd_calls_original_set = set(bnd_calls_original)
						bnd_intersection = bnd_calls_original_set.intersection(bnd_calls_shuffled_set)
						bnd_calls_original_only = bnd_calls_original_set.difference(bnd_calls_shuffled_set)
						bnd_calls_shuffled_only = bnd_calls_shuffled_set.difference(bnd_calls_original_set)

						bedfile_intersection_statistics = [shuffled_filter, variant_type, len(bnd_intersection), len(bnd_calls_original_only ), len(bnd_calls_shuffled_only)]
						intersection_statistics.append(bedfile_intersection_statistics)

					if len(bnd_calls_original_only) > 0:
						for x in bnd_calls_original_only:
							x_split = x.split("\t")
							unique_entry = ["Original VCF", shuffled_filter, variant_type, x_split[0], x_split[1], x_split[2], "NA"]
							unique_calls_original.append(','.join(str(e) for e in unique_entry))

					if bnd_calls_original_only is not None	and	len(bnd_calls_shuffled_only) > 0:
						for x in bnd_calls_shuffled_only:
							x_split = x.split("\t")
							unique_entry = ["Shuffled VCF", shuffled_filter, variant_type, x_split[0], x_split[1], x_split[2], "NA"]
							unique_calls_shuffled.append(','.join(str(e) for e in unique_entry))

	else:
		for shuffled_filter in svim_shuffled.keys():
			shuffled_variant_types = svim_shuffled[shuffled_filter].keys()
			original_variant_types = svim_original[shuffled_filter].keys()
			unique_variant_types = set(shuffled_variant_types) ^ set(original_variant_types) # Identify any variant types only found in the original or shuffled vcf files

			if len(unique_variant_types) != 0:
				for svtype in unique_variant_types:
					if svtype not in original_variant_types:
						warning_messages.add(("Variant type: " + svtype + " not found in the original vcf file for the filter type: " + shuffled_filter))
					if svtype not in shuffled_variant_types:
						warning_messages.add(("Variant type: " + svtype + " not found in the shuffled vcf file for the filter type: " + shuffled_filter))

			for variant_type in shuffled_variant_types:
				if variant_type != "BND":
					sv_calls_shuffled = svim_shuffled[shuffled_filter][variant_type]
					bedlines_shuffled = '\n'.join(str(e) for e in sv_calls_shuffled)
					bedfile_shuffled = BedTool(bedlines_shuffled, from_string=True)

					if variant_type not in original_variant_types:
						shuffled_only = bedfile_shuffled
						bedfile_intersection_statistics = [shuffled_filter, variant_type, 0, 0, len(shuffled_only)]
						original_only = None
						#print(bedfile_intersection_statistics)
						intersection_statistics.append(bedfile_intersection_statistics)
					else:

						sv_calls_original = svim_original[shuffled_filter][variant_type]
						bedlines_original = '\n'.join(str(e) for e in sv_calls_original)
						bedfile_original = BedTool(bedlines_original, from_string=True)
						bedfile_intersection = bedfile_original.intersect(bedfile_shuffled)
						original_only = bedfile_original.intersect(bedfile_shuffled, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
						shuffled_only = bedfile_shuffled.intersect(bedfile_original, v=True) # -v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
						bedfile_intersection_statistics = [shuffled_filter, variant_type, len(bedfile_intersection), len(original_only), len(shuffled_only)]
						#print(bedfile_intersection_statistics)
						intersection_statistics.append(bedfile_intersection_statistics)

					if original_only is not None and len(original_only) > 0:
						for x in original_only:
							unique_entry = ["Original VCF", shuffled_filter, variant_type, x.chrom, x.start, x.stop, len(x)]
							unique_calls_original.append(','.join(str(e) for e in unique_entry))

					if len(shuffled_only) > 0:
						for x in shuffled_only:
							unique_entry = ["Shuffled VCF", shuffled_filter, variant_type, x.chrom, x.start, x.stop, len(x)]
							unique_calls_shuffled.append(','.join(str(e) for e in unique_entry))
				else:
					bnd_calls_shuffled = svim_shuffled[shuffled_filter][variant_type]
					bnd_calls_shuffled_set = set(bnd_calls_shuffled)

					if variant_type not in original_variant_types:
						bnd_calls_shuffled_only = bnd_calls_shuffled
						bedfile_intersection_statistics = [shuffled_filter, variant_type, 0, 0, len(bnd_calls_shuffled_only)]
						bnd_calls_original_only = None
						#print(bedfile_intersection_statistics)
						intersection_statistics.append(bedfile_intersection_statistics)

					else:
						bnd_calls_original = svim_original[shuffled_filter][variant_type]
						bnd_calls_original_set = set(bnd_calls_original)
						bnd_intersection = bnd_calls_original_set.intersection(bnd_calls_shuffled_set)
						bnd_calls_original_only = bnd_calls_original_set.difference(bnd_calls_shuffled_set)
						bnd_calls_shuffled_only = bnd_calls_shuffled_set.difference(bnd_calls_original_set)

						bedfile_intersection_statistics = [shuffled_filter, variant_type, len(bnd_intersection), len(bnd_calls_original_only ), len(bnd_calls_shuffled_only)]
						intersection_statistics.append(bedfile_intersection_statistics)

					if bnd_calls_original_only is not None and len(bnd_calls_original_only) > 0:
						for x in bnd_calls_original_only:
							x_split = x.split("\t")
							unique_entry = ["Original VCF", shuffled_filter, variant_type, x_split[0], x_split[1], x_split[2], "NA"]
							unique_calls_original.append(','.join(str(e) for e in unique_entry))

					if len(bnd_calls_shuffled_only) > 0:
						for x in bnd_calls_shuffled_only:
							x_split = x.split("\t")
							unique_entry = ["Shuffled VCF", shuffled_filter, variant_type, x_split[0], x_split[1], x_split[2], "NA"]
							unique_calls_shuffled.append(','.join(str(e) for e in unique_entry))

	#for row in intersection_statistics[1:]:
	#	print(row)
	#	print(int(row[3]))
	for x in warning_messages:
		print(x)

	total_intersection = sum([int(row[2]) for row in intersection_statistics[1:]])
	total_original_only = sum([int(row[3]) for row in intersection_statistics[1:]])
	#print(total_original_only)
	total_shuffled_only = sum([int(row[4]) for row in intersection_statistics[1:]])
	intersection_statistics.append(["Total","ALL_SVS", total_intersection, total_original_only, total_shuffled_only])
	#for x in intersection_statistics:
	#	print(x)

	return intersection_statistics, unique_calls_original, unique_calls_shuffled

def compare_vcf_breakpoints(svim_shuffled, svim_original):
	intersection_statistics = [["FILTER,SVTYPE,INTERSECTION,SAME_BREAKPOINTS,DIFFERENT_BREAKPOINTS"]] # Store comparison statistics: shuffled_filter, variant_type, intersection size, unique to original, unique to shuffled
	different_breakpoints = [] # Store variants called by both the original and shuffled VCF files but with different breakpoints
	shuffled_filters = svim_shuffled.keys()
	original_filters = svim_original.keys()
	unique_filters = set(shuffled_filters) ^ set(original_filters)

	if len(unique_filters) == 0:
		for shuffled_filter in svim_shuffled.keys():
			shuffled_variant_types = svim_shuffled[shuffled_filter].keys()
			original_variant_types = svim_original[shuffled_filter].keys()
			unique_variant_types = set(shuffled_variant_types) ^ set(original_variant_types) # Identify any variant types only found in the original or shuffled vcf files

			if len(unique_variant_types) == 0:
				for variant_type in shuffled_variant_types:
					if variant_type != "BND" and variant_type != "INS":
						same_breakpoint_count = 0
						different_breakpoint_count = 0
						sv_calls_shuffled = svim_shuffled[shuffled_filter][variant_type]
						bedlines_shuffled = '\n'.join(str(e) for e in sv_calls_shuffled)
						bedfile_shuffled = BedTool(bedlines_shuffled, from_string=True)
						sv_calls_original = svim_original[shuffled_filter][variant_type]
						bedlines_original = '\n'.join(str(e) for e in sv_calls_original)
						bedfile_original = BedTool(bedlines_original, from_string=True)
						bedfile_intersection = bedfile_original.intersect(bedfile_shuffled, wa=True, wb=True)

						for x in bedfile_intersection:
							chromosome = x[0]
							original_start = x[1]
							original_end = x[2]
							shuffled_start = x[4]
							shuffled_end = x[5]
							if original_start != shuffled_start or original_end != shuffled_end:
								different_breakpoint_count += 1
							else:
								same_breakpoint_count +=1
						intersection_statistics_line = [shuffled_filter,variant_type, len(bedfile_intersection),same_breakpoint_count,different_breakpoint_count]
						intersection_statistics.append(intersection_statistics_line)
						#intersection_statistics.append(','.join(str(e) for e in intersection_statistics_line))
	total_intersection = sum([int(row[2]) for row in intersection_statistics[1:]])
	total_same_breakpoints = sum([int(row[3]) for row in intersection_statistics[1:]])
	total_different_breakpoints = sum([int(row[4]) for row in intersection_statistics[1:]])
	intersection_statistics.append(["Total","ALL_SVS", total_intersection, total_same_breakpoints, total_different_breakpoints])

	return intersection_statistics, different_breakpoints

# Write variant calls to disk
def write_results(intersection_statistics_total_calls, unique_calls_original, unique_calls_shuffled,intersection_statistics_breakpoints, different_breakpoints):
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
	# Check if files exist. If not, write to file.
	if os.path.isfile(output_file):
		print(str(output_file) + " exists. Not overwriting.")
	else:
		intersection_statistics_total_calls_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
		for line in intersection_statistics_total_calls:
			intersection_statistics_total_calls_string.append(','.join(str(e) for e in line))
		with open(output_file, 'w', newline='\n') as f:
			f.writelines("%s\n" % l for l in intersection_statistics_total_calls_string)

	output_file = outdir + "/summary_breakpoints.csv"
	# Check if files exist. If not, write to file.
	if os.path.isfile(output_file):
		print(str(output_file) + " exists. Not overwriting.")
	else:
		intersection_statistics_breakpoints_string = [] # Stats were stored in a list containing a list for each row. Convert the rows to strings.
		for line in intersection_statistics_breakpoints:
			intersection_statistics_breakpoints_string.append(','.join(str(e) for e in line))
		with open(output_file, 'w', newline='\n') as f:
			f.writelines("%s\n" % l for l in intersection_statistics_breakpoints_string)

if args.sv_caller.lower() == "svim":
	with open(args.shuffled_vcf) as f:
		variants = f.readlines()
	svim_shuffled = parse_svim(variants)

	with open(args.original_vcf) as f:
		variants = f.readlines()
	svim_original = parse_svim(variants)
	vcf_comparison_total_svs = compare_vcf_total_predictions(svim_shuffled, svim_original)
	vcf_comparison_breakpoints = compare_vcf_breakpoints(svim_shuffled, svim_original)
	write_results(vcf_comparison_total_svs[0], vcf_comparison_total_svs[1], vcf_comparison_total_svs[2],vcf_comparison_breakpoints[0], vcf_comparison_breakpoints[1])

if args.sv_caller.lower() == "sniffles":
	with open(args.shuffled_vcf) as f:
		variants = f.readlines()
	sniffles_shuffled = parse_sniffles(variants)

	with open(args.original_vcf) as f:
		variants = f.readlines()
	sniffles_original = parse_sniffles(variants)
	vcf_comparison_total_svs = compare_vcf_total_predictions(sniffles_shuffled, sniffles_original)
	vcf_comparison_breakpoints = compare_vcf_breakpoints(sniffles_shuffled, sniffles_original)
	write_results(vcf_comparison_total_svs[0], vcf_comparison_total_svs[1], vcf_comparison_total_svs[2],vcf_comparison_breakpoints[0], vcf_comparison_breakpoints[1])

if args.sv_caller.lower() == "pbsv":
	with open(args.shuffled_vcf) as f:
		variants = f.readlines()
	pbsv_shuffled = parse_pbsv(variants)

	with open(args.original_vcf) as f:
		variants = f.readlines()
	pbsv_original = parse_pbsv(variants)
	vcf_comparison_total_svs = compare_vcf_total_predictions(pbsv_shuffled, pbsv_original)
	vcf_comparison_breakpoints = compare_vcf_breakpoints(pbsv_shuffled, pbsv_original)
	write_results(vcf_comparison_total_svs[0], vcf_comparison_total_svs[1], vcf_comparison_total_svs[2],vcf_comparison_breakpoints[0], vcf_comparison_breakpoints[1])
