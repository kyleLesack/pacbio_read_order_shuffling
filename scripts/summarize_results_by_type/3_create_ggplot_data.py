#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 20:45:14 2022

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
INTERSECTION_FULL_DEPTH_FILE = "4_results/full_depth/sv_intersection_agreement/CALLER/CALLER-ALIGNER.csv"
BREAKPOINT_SUBSAMPLED_FILE = "4_results/subsampled/DEPTH/breakpoint_agreement/CALLER/CALLER-ALIGNER.csv"
INTERSECTION_SUBSAMPLED_FILE = "4_results/subsampled/DEPTH/sv_intersection_agreement/CALLER/CALLER-ALIGNER.csv"

def add_totals_breakpoint(input_file, caller, aligner, depth):
	df=pd.read_csv(input_file)
	sum_row = {col: df[col].sum() for col in df}
	sum_df = pd.DataFrame(sum_row, index=["Total"])
	sum_df['SV TYPE'] = "ALL SVs"
	df = pd.concat([df,sum_df])
	df['Discordant (%)'] = df['Different breakpoints'] / df['Intersection']
	df['Discordant (%)'] = df['Discordant (%)'] * 100
	df2 = df.loc[:,["Different breakpoints","Same breakpoints", "Discordant (%)"]]
	stacked = df2.stack()
	stacked_df = pd.Series.to_frame(stacked)
	stacked_df.index = stacked_df.index.set_names(['SV Type', 'Agreement'])
	stacked_df.columns =['Count']
	stacked_df.reset_index(inplace=True)
	stacked_df['Caller'] = caller
	stacked_df['Aligner'] = aligner
	if depth is not None:
		stacked_df['Depth'] = depth
	total_df = stacked_df.iloc[-3:]
	return(total_df)

def add_totals_Intersection(input_file, caller, aligner, depth):
	df=pd.read_csv(input_file, index_col="SV TYPE")
	sum_row = {col: df[col].sum() for col in df}
	sum_df = pd.DataFrame(sum_row, index=["Total"])
	df = pd.concat([df,sum_df])
	df['Unique'] = df['Original only'] + df['Shuffled only']
	df['Total_SVs'] = df['Unique'] + df['Intersection']
	df['Non-Overlapping (%)'] = df['Unique'] / df['Total_SVs']
	df['Non-Overlapping (%)'] = df['Non-Overlapping (%)'] * 100
	df2 = df.loc[:,["Intersection","Unique","Total_SVs","Non-Overlapping (%)"]]
	stacked = df2.stack()
	stacked_df = pd.Series.to_frame(stacked)
	stacked_df.index = stacked_df.index.set_names(['SV TYPE', 'Agreement'])
	stacked_df.columns =['Count']
	stacked_df.reset_index(inplace=True)
	stacked_df['Caller'] = caller
	stacked_df['Aligner'] = aligner
	if depth is not None:
		stacked_df['Depth'] = depth
	total_df = stacked_df.iloc[-4:]
	return(total_df)

# Write results to disk
def write_results(result_df, outfile):
	parent_path = outfile.rpartition('/')[0]
	if not os.path.exists(parent_path):
		os.makedirs(parent_path)

	result_df.to_csv(outfile + ".tsv", sep='\t', index = False, float_format="%.0f")
	result_df.to_csv(outfile + ".csv", sep=',', index = False, float_format="%.0f")

# Add totals and percentage Discordant to breakpoint results
for caller in CALLERS:
	if caller == "pbsv":
		for aligner in PBSV_ALIGNERS:
			breakpoint_file = BREAKPOINT_FULL_DEPTH_FILE
			breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_breakpoint(breakpoint_file, caller, aligner, None)
			outpath = breakpoint_file.rpartition('/')[0] + "/ggplot/"
			outfile = breakpoint_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

	if caller == "sniffles":
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			breakpoint_file = BREAKPOINT_FULL_DEPTH_FILE
			breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_breakpoint(breakpoint_file, caller, aligner, None)
			outpath = breakpoint_file.rpartition('/')[0] + "/ggplot/"
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
			result_df = add_totals_breakpoint(breakpoint_file, "svim", aligner, None)
			outpath = breakpoint_file.rpartition('/')[0] + "/ggplot/"
			outfile = breakpoint_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

# Add totals to Intersection results
for caller in CALLERS:
	if caller == "pbsv":
		for aligner in PBSV_ALIGNERS:
			sv_Intersection_file = INTERSECTION_FULL_DEPTH_FILE
			sv_Intersection_file = sv_Intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_Intersection(sv_Intersection_file, caller, aligner, None)
			outpath = sv_Intersection_file.rpartition('/')[0] + "/ggplot/"
			outfile = sv_Intersection_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

	if caller == "sniffles":
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			sv_Intersection_file = INTERSECTION_FULL_DEPTH_FILE
			sv_Intersection_file = sv_Intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			result_df = add_totals_Intersection(sv_Intersection_file, caller, aligner, None)
			outpath = sv_Intersection_file.rpartition('/')[0] + "/ggplot/"
			outfile = sv_Intersection_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

	if caller == "svim/qual_QUAL/":
		#for qual in ["0","15"]:
		qual = "15"
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			sv_Intersection_file = INTERSECTION_FULL_DEPTH_FILE
			sv_Intersection_file = sv_Intersection_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("QUAL", qual).replace("CALLER", "svim")
			result_df = add_totals_Intersection(sv_Intersection_file, "svim", aligner, None)
			outpath = sv_Intersection_file.rpartition('/')[0] + "/ggplot/"
			outfile = sv_Intersection_file.rpartition('/')[2]
			outfile = outpath + outfile
			outfile = outfile.rpartition('.')[0]
			write_results(result_df, outfile)

# Subsampled analysis
# Add totals and percentage Discordant to breakpoint results
for caller in CALLERS:
	for depth in SUBSAMPLE_DEPTHS:
		if caller == "pbsv":
			for aligner in PBSV_ALIGNERS:
				breakpoint_file = BREAKPOINT_SUBSAMPLED_FILE
				breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_breakpoint(breakpoint_file, caller, aligner, depth)
				outpath = breakpoint_file.rpartition('/')[0] + "/ggplot/"
				outfile = breakpoint_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

		if caller == "sniffles":
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				breakpoint_file = BREAKPOINT_SUBSAMPLED_FILE
				breakpoint_file = breakpoint_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_breakpoint(breakpoint_file, caller, aligner, depth)
				outpath = breakpoint_file.rpartition('/')[0] + "/ggplot/"
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
				result_df = add_totals_breakpoint(breakpoint_file, "svim", aligner, depth)
				outpath = breakpoint_file.rpartition('/')[0] + "/ggplot/"
				outfile = breakpoint_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

# Add totals to Intersection results
for caller in CALLERS:
	for depth in SUBSAMPLE_DEPTHS:
		if caller == "pbsv":
			for aligner in PBSV_ALIGNERS:
				sv_Intersection_file = INTERSECTION_SUBSAMPLED_FILE
				sv_Intersection_file = sv_Intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_Intersection(sv_Intersection_file, caller, aligner, depth)
				outpath = sv_Intersection_file.rpartition('/')[0] + "/ggplot/"
				outfile = sv_Intersection_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

		if caller == "sniffles":
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				sv_Intersection_file = INTERSECTION_SUBSAMPLED_FILE
				sv_Intersection_file = sv_Intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				result_df = add_totals_Intersection(sv_Intersection_file, caller, aligner, depth)
				outpath = sv_Intersection_file.rpartition('/')[0] + "/ggplot/"
				outfile = sv_Intersection_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)

		if caller == "svim/qual_QUAL/":
			#for qual in ["0","15"]:
			qual = "15"
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				sv_Intersection_file = INTERSECTION_SUBSAMPLED_FILE
				sv_Intersection_file = sv_Intersection_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("QUAL", qual).replace("CALLER", "svim").replace("DEPTH", depth)
				result_df = add_totals_Intersection(sv_Intersection_file, "svim", aligner, depth)
				outpath = sv_Intersection_file.rpartition('/')[0] + "/ggplot/"
				outfile = sv_Intersection_file.rpartition('/')[2]
				outfile = outpath + outfile
				outfile = outfile.rpartition('.')[0]
				write_results(result_df, outfile)
