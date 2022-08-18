#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 22:06:48 2022

@author: kyle
"""

import pandas as pd
import os
#os.chdir("/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/") # Change working directory to parent dir to the results

SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["pbmm2"]
PBSV_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV"]
SNIFFLES_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV"]
SVIM_SV_TYPES = ["BND", "DEL", "DUP:INT","DUP:TANDEM", "INS","INV"]
SUBSAMPLE_DEPTHS = ["10X","20X","40X","60X"]
CALLERS = ["pbsv","sniffles","svim/qual_QUAL/"]

BREAKPOINT_FULL_DEPTH_FILE = "4_results/full_depth/breakpoint_agreement/CALLER/CALLER-ALIGNER.csv"
Intersection_FULL_DEPTH_FILE = "4_results/full_depth/sv_intersection_agreement/CALLER/CALLER-ALIGNER.csv"
BREAKPOINT_SUBSAMPLED_FILE = "4_results/subsampled/DEPTH/breakpoint_agreement/CALLER/CALLER-ALIGNER.csv"
Intersection_SUBSAMPLED_FILE = "4_results/subsampled/DEPTH/sv_intersection_agreement/CALLER/CALLER-ALIGNER.csv"

def add_totals_breakpoint(input_file):
	df=pd.read_csv(input_file)
	sum_row = {col: df[col].sum() for col in df}
	sum_df = pd.DataFrame(sum_row, index=["Total"])
	sum_df['SV TYPE'] = "ALL SVs"
	df = pd.concat([df,sum_df])
	df['Discordant (%)'] = df['Different breakpoints'] / df['Intersection']
	df['Discordant (%)'] = df['Discordant (%)'] * 100
	return(df)

def add_totals_intersection(input_file):
	df=pd.read_csv(input_file)
	sum_row = {col: df[col].sum() for col in df}
	sum_df = pd.DataFrame(sum_row, index=["Total"])
	sum_df['SV TYPE'] = "ALL SVs"
	df = pd.concat([df,sum_df])
	df['Unique'] = df['Original only'] + df['Shuffled only']
	return(df)

# Write results to disk
def write_results(result_df, outfile):
	parent_path = outfile.rpartition('/')[0]
	if not os.path.exists(parent_path):
		os.makedirs(parent_path)

	result_df.to_csv(outfile + ".tsv", sep='\t', index = False, float_format="%.0f")
	result_df.to_csv(outfile + ".csv", sep=',', index = False, float_format="%.0f")

# Add totals and percentage Discordant (%) to breakpoint results
for caller in CALLERS:
	if caller == "pbsv":
		for aligner in PBSV_ALIGNERS:
			breakpoint_file = BREAKPOINT_FULL_DEPTH_FILE
			breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_breakpoint(breakpoint_file)
			outpath = breakpoint_file.rpartition('/')[0] + "/totals/"
			outfile = breakpoint_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

	if caller == "sniffles":
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			breakpoint_file = BREAKPOINT_FULL_DEPTH_FILE
			breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_breakpoint(breakpoint_file)
			outpath = breakpoint_file.rpartition('/')[0] + "/totals/"
			outfile = breakpoint_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

	if caller == "svim/qual_QUAL/":
		#for qual in ["0","15"]:
		qual = "15"
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			breakpoint_file = BREAKPOINT_FULL_DEPTH_FILE
			breakpoint_file = breakpoint_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("QUAL", qual).replace("CALLER", "svim")
			result_df = add_totals_breakpoint(breakpoint_file)
			outpath = breakpoint_file.rpartition('/')[0] + "/totals/"
			outfile = breakpoint_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

# Add totals to intersection results
for caller in CALLERS:
	if caller == "pbsv":
		for aligner in PBSV_ALIGNERS:
			sv_intersection_file = Intersection_FULL_DEPTH_FILE
			sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_intersection(sv_intersection_file)
			outpath = sv_intersection_file.rpartition('/')[0] + "/totals/"
			outfile = sv_intersection_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

	if caller == "sniffles":
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			sv_intersection_file = Intersection_FULL_DEPTH_FILE
			sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_intersection(sv_intersection_file)
			outpath = sv_intersection_file.rpartition('/')[0] + "/totals/"
			outfile = sv_intersection_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

	if caller == "svim/qual_QUAL/":
		#for qual in ["0","15"]:
		qual = "15"
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			sv_intersection_file = Intersection_FULL_DEPTH_FILE
			sv_intersection_file = sv_intersection_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("QUAL", qual).replace("CALLER", "svim")
			result_df = add_totals_intersection(sv_intersection_file)
			outpath = sv_intersection_file.rpartition('/')[0] + "/totals/"
			outfile = sv_intersection_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

# Subsampled analysis
# Add totals and percentage Discordant (%) to breakpoint results
for caller in CALLERS:
	for depth in SUBSAMPLE_DEPTHS:
		if caller == "pbsv":
			for aligner in PBSV_ALIGNERS:
				breakpoint_file = BREAKPOINT_SUBSAMPLED_FILE
				breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_breakpoint(breakpoint_file)
				outpath = breakpoint_file.rpartition('/')[0] + "/totals/"
				outfile = breakpoint_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

		if caller == "sniffles":
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				breakpoint_file = BREAKPOINT_SUBSAMPLED_FILE
				breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_breakpoint(breakpoint_file)
				outpath = breakpoint_file.rpartition('/')[0] + "/totals/"
				outfile = breakpoint_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

		if caller == "svim/qual_QUAL/":
			#for qual in ["0","15"]:
			qual = "15"
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				breakpoint_file = BREAKPOINT_SUBSAMPLED_FILE
				breakpoint_file = breakpoint_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("QUAL", qual).replace("CALLER", "svim").replace("DEPTH", depth)
				result_df = add_totals_breakpoint(breakpoint_file)
				outpath = breakpoint_file.rpartition('/')[0] + "/totals/"
				outfile = breakpoint_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

# Add totals to intersection results for subsampled data
for caller in CALLERS:
	for depth in SUBSAMPLE_DEPTHS:
		if caller == "pbsv":
			for aligner in PBSV_ALIGNERS:
				sv_intersection_file = Intersection_SUBSAMPLED_FILE
				sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_intersection(sv_intersection_file)
				outpath = sv_intersection_file.rpartition('/')[0] + "/totals/"
				outfile = sv_intersection_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]

				write_results(result_df, outfile)

		if caller == "sniffles":
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				sv_intersection_file = Intersection_SUBSAMPLED_FILE
				sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_intersection(sv_intersection_file)
				outpath = sv_intersection_file.rpartition('/')[0] + "/totals/"
				outfile = sv_intersection_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

		if caller == "svim/qual_QUAL/":
			#for qual in ["0","15"]:
			qual = "15"
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				sv_intersection_file = Intersection_SUBSAMPLED_FILE
				sv_intersection_file = sv_intersection_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("QUAL", qual).replace("CALLER", "svim").replace("DEPTH", depth)
				result_df = add_totals_intersection(sv_intersection_file)
				outpath = sv_intersection_file.rpartition('/')[0] + "/totals/"
				outfile = sv_intersection_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)
