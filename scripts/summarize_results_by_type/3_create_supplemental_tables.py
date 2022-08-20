#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 15:57:39 2022

@author: kyle
"""
import pandas as pd
import os
from collections import defaultdict
from pathlib import Path  

os.chdir("/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling/") # Change working directory to parent dir to the results

SNIFFLES_SVIM_ALIGNERS = ["minimap2", "ngmlr", "pbmm2"]
PBSV_ALIGNERS = ["pbmm2"]
PBSV_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV"]
SNIFFLES_SV_TYPES = ["BND", "DEL", "DUP", "INS","INV"]
SVIM_SV_TYPES = ["BND", "DEL", "DUP:INT","DUP:TANDEM", "INS","INV"]
SUBSAMPLE_DEPTHS = ["10X","20X","40X","60X"]
CALLERS = ["pbsv","sniffles","svim/qual_15/"]

BREAKPOINT_FULL_DEPTH_FILE = "4_results/full_depth/breakpoint_agreement/CALLER/CALLER-ALIGNER.csv"
INTERSECTION_FULL_DEPTH_FILE = "4_results/full_depth/sv_intersection_agreement/CALLER/CALLER-ALIGNER.csv"
BREAKPOINT_SUBSAMPLED_FILE = "4_results/subsampled/DEPTH/breakpoint_agreement/CALLER/CALLER-ALIGNER.csv"
INTERSECTION_SUBSAMPLED_FILE = "4_results/subsampled/DEPTH/sv_intersection_agreement/CALLER/CALLER-ALIGNER.csv"
EXCEL_FILEPATH = Path('4_results/full_depth/sv_intersection_agreement/supplemental/sv_agreement.xlsx')

def write_excel(filepath,df,sheetname):
	if filepath.is_file():
		print("Appending to: " + str(filepath))
		with pd.ExcelWriter(filepath, mode='a') as writer:  
			df.to_excel(writer, sheet_name=sheetname)
	else:
		print("Writing to: " + str(filepath))
		with pd.ExcelWriter(filepath, mode='w') as writer:  
			df.to_excel(writer, sheet_name=sheetname)	

def merge_sv_agreement_dfs(caller_aligner_sv_agreement_dict, caller):
	combined_df = pd.DataFrame()
	for aligner in caller_aligner_sv_agreement_dict.keys():
		for df in caller_aligner_sv_agreement_dict[aligner]:
			df['Aligner'] = aligner
			combined_df = pd.concat([combined_df, df], axis=0)
	
	intersection_df = combined_df[["SV TYPE", "Aligner", "Intersection"]]
	intersection_df = intersection_df.reset_index(drop=True)
	intersection_df_pivot = intersection_df.pivot(index="SV TYPE", columns="Aligner", values="Intersection")
	intersection_df_pivot = intersection_df_pivot.add_prefix(caller)
	intersection_df_pivot = intersection_df_pivot.add_suffix('_int')
	
	unique_df = combined_df[["SV TYPE", "Aligner", "Unique"]]
	unique_df = unique_df.reset_index(drop=True)
	unique_df_pivot = unique_df.pivot(index="SV TYPE", columns="Aligner", values="Unique")
	unique_df_pivot = unique_df_pivot.add_prefix(caller)
	unique_df_pivot = unique_df_pivot.add_suffix('_uni')

	unique_proportion_df = combined_df[["SV TYPE", "Aligner", "Unique proportion"]]
	unique_proportion_df = unique_proportion_df.reset_index(drop=True)
	unique_proportion_df_pivot = unique_proportion_df.pivot(index="SV TYPE", columns="Aligner", values="Unique proportion")
	unique_proportion_df_pivot = unique_proportion_df_pivot.add_prefix(caller)
	unique_proportion_df_pivot = unique_proportion_df_pivot.add_suffix('_uniq_prop')
	intersection_unique_df = pd.concat([intersection_df_pivot, unique_df_pivot], axis=1)
	intersection_unique_uniq_prop_df = pd.concat([intersection_unique_df, unique_proportion_df_pivot], axis=1)
	#intersection_unique_df = intersection_unique_df.rename_axis(index=("SV Type"), columns=None)
	
	idx = intersection_unique_df.index.tolist()
	idx.append(idx.pop(0)) # Move first element to last
	intersection_unique_uniq_prop_df = intersection_unique_uniq_prop_df.reindex(idx)
	unique_proportion_df_pivot = unique_proportion_df_pivot.reindex(idx)
	
	#intersection_unique_df = intersection_unique_df.reindex(idx)
	#intersection_unique_df = intersection_unique_df.reindex(sorted(intersection_unique_df.columns), axis=1) # Sort columns so intersection and unique values columns are adjacent for each aligner	
	intersection_unique_uniq_prop_df = intersection_unique_uniq_prop_df.reindex(sorted(intersection_unique_uniq_prop_df.columns), axis=1) # Sort columns so intersection and unique values columns are adjacent for each aligner	
	
	return(intersection_unique_uniq_prop_df,unique_proportion_df_pivot)

def merge_sv_agreement_proportions(sv_agreement_uniq_prop_df_list):
	combined_df = pd.DataFrame()
	for x in sv_agreement_uniq_prop_df_list:
		df = x[1]
		df.rename(index={'DUP:TANDEM': 'DUP'},inplace=True)
		combined_df = pd.concat([combined_df, df], axis=1)
	idx = combined_df.index.tolist()
	idx.sort() 
	idx.append(idx.pop(0)) # Move first element to last
	combined_df = combined_df.reindex(idx)

	return(combined_df)

# If excel files exist, remove them. If they already exist, appending to a file with same sheet names causes an error.
filepath = EXCEL_FILEPATH
if filepath.is_file():
	os.remove(filepath)

sv_agreement_uniq_prop_df_list = [] # List to store dataframes with the unique call proportions for each caller/aligner combo

# Open SV intersection files
for caller in CALLERS:
	
	if caller == "pbsv":
		pbsv_aligner_sv_agreement_dict = defaultdict(list)
		for aligner in PBSV_ALIGNERS:
			sv_intersection_file = INTERSECTION_FULL_DEPTH_FILE
			sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			sv_intersection_file_df = pd.read_csv(sv_intersection_file)
			pbsv_aligner_sv_agreement_dict[aligner].append(sv_intersection_file_df)
		pbsv_merged = merge_sv_agreement_dfs(pbsv_aligner_sv_agreement_dict,"pb_")
		sv_agreement_uniq_prop_df_list.append(("pbsv",pbsv_merged[1]))
		filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S1-sv_agreement_pbsv.csv')  
		filepath.parent.mkdir(parents=True, exist_ok=True)  
		pbsv_merged[0].to_csv(filepath, float_format='%.2f')  
		write_excel(EXCEL_FILEPATH,pbsv_merged[0],"S1_pbsv")

	if caller == "sniffles":
		#sniffles_aligner_sv_agreement = []
		sniffles_aligner_sv_agreement_dict = defaultdict(list)
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			sv_intersection_file = INTERSECTION_FULL_DEPTH_FILE
			sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner)
			sv_intersection_file_df = pd.read_csv(sv_intersection_file)
			sniffles_aligner_sv_agreement_dict[aligner].append(sv_intersection_file_df)
		sniffles_merged = merge_sv_agreement_dfs(sniffles_aligner_sv_agreement_dict,"sn_")
		sv_agreement_uniq_prop_df_list.append(("sniffles",sniffles_merged[1]))
		filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S2-sv_agreement_sniffles.csv')  
		filepath.parent.mkdir(parents=True, exist_ok=True)  
		sniffles_merged[0].to_csv(filepath, float_format='%.2f')  
		write_excel(EXCEL_FILEPATH,sniffles_merged[0],"S2_Sniffles")
					
	if caller == "svim/qual_15/":
		svim_aligner_sv_agreement_dict = defaultdict(list)
		for aligner in SNIFFLES_SVIM_ALIGNERS:
			sv_intersection_file = INTERSECTION_FULL_DEPTH_FILE
			sv_intersection_file = sv_intersection_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("CALLER", "svim")
			sv_intersection_file_df = pd.read_csv(sv_intersection_file)
			svim_aligner_sv_agreement_dict[aligner].append(sv_intersection_file_df)
		svim_merged = merge_sv_agreement_dfs(svim_aligner_sv_agreement_dict,"sv_")
		sv_agreement_uniq_prop_df_list.append(("svim",svim_merged[1]))
		filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S3-sv_agreement_svim.csv') 
		filepath.parent.mkdir(parents=True, exist_ok=True)  
		svim_merged[0].to_csv(filepath, float_format='%.2f')  
		write_excel(EXCEL_FILEPATH,svim_merged[0],"S3_SVIM")

uniq_prop_all_callers = merge_sv_agreement_proportions(sv_agreement_uniq_prop_df_list)
filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S4-uniq_prop_all_callers.csv') 
uniq_prop_all_callers.to_csv(filepath, float_format='%.2f') 
write_excel(EXCEL_FILEPATH,uniq_prop_all_callers,"S4")

# Subsampled analysis


## Open SV intersection files

sv_agreement_depth_list = []
for depth in SUBSAMPLE_DEPTHS:	
	sv_agreement_caller_list = [] # List to store dataframes with the unique call proportions for each caller/aligner combo
	for caller in CALLERS:	
		if caller == "pbsv":
			pbsv_aligner_sv_agreement_dict = defaultdict(list)
			for aligner in PBSV_ALIGNERS:
				sv_intersection_file = INTERSECTION_SUBSAMPLED_FILE
				sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				sv_intersection_file_df = pd.read_csv(sv_intersection_file)
				pbsv_aligner_sv_agreement_dict[aligner].append(sv_intersection_file_df)
			pbsv_merged = merge_sv_agreement_dfs(pbsv_aligner_sv_agreement_dict,"pbsv_")
			sv_agreement_caller_list.append(("pbsv",pbsv_merged[1]))
			#filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S1-sv_agreement_pbsv.csv')  
			#filepath.parent.mkdir(parents=True, exist_ok=True)  
			#pbsv_merged[0].to_csv(filepath, float_format='%.2f')  
			#write_excel(EXCEL_FILEPATH,pbsv_merged[0],"S1_pbsv")
	
		if caller == "sniffles":
			#sniffles_aligner_sv_agreement = []
			sniffles_aligner_sv_agreement_dict = defaultdict(list)
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				sv_intersection_file = INTERSECTION_SUBSAMPLED_FILE
				sv_intersection_file = sv_intersection_file.replace("CALLER", caller).replace("ALIGNER", aligner).replace("DEPTH", depth)
				sv_intersection_file_df = pd.read_csv(sv_intersection_file)
				sniffles_aligner_sv_agreement_dict[aligner].append(sv_intersection_file_df)
			sniffles_merged = merge_sv_agreement_dfs(sniffles_aligner_sv_agreement_dict,"Sniffles_")
			sv_agreement_caller_list.append(("sniffles",sniffles_merged[1]))
			#filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S2-sv_agreement_sniffles.csv')  
			#filepath.parent.mkdir(parents=True, exist_ok=True)  
			#sniffles_merged[0].to_csv(filepath, float_format='%.2f')  
			#write_excel(EXCEL_FILEPATH,sniffles_merged[0],"S2_Sniffles")
						
		if caller == "svim/qual_15/":
			svim_aligner_sv_agreement_dict = defaultdict(list)
			for aligner in SNIFFLES_SVIM_ALIGNERS:
				sv_intersection_file = INTERSECTION_SUBSAMPLED_FILE
				sv_intersection_file = sv_intersection_file.replace("CALLER", caller,1).replace("ALIGNER", aligner).replace("CALLER", "svim").replace("DEPTH", depth)
				sv_intersection_file_df = pd.read_csv(sv_intersection_file)
				svim_aligner_sv_agreement_dict[aligner].append(sv_intersection_file_df)
			svim_merged = merge_sv_agreement_dfs(svim_aligner_sv_agreement_dict,"SVIM_")
			sv_agreement_caller_list.append(("svim",svim_merged[1]))
			#filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S3-sv_agreement_svim.csv') 
			#filepath.parent.mkdir(parents=True, exist_ok=True)  
			#svim_merged[0].to_csv(filepath, float_format='%.2f')  
			#write_excel(EXCEL_FILEPATH,svim_merged[0],"S3_SVIM")

	uniq_prop_all_callers = merge_sv_agreement_proportions(sv_agreement_caller_list)
	#uniq_prop_all_callers['Depth'] = depth
	depth_column = [depth] *len(uniq_prop_all_callers)
	uniq_prop_all_callers.insert(0, 'Depth', depth_column)
	
	sv_agreement_depth_list.append(uniq_prop_all_callers)

combined_df = pd.DataFrame()
for df in sv_agreement_depth_list:
	#print(x[1])
	#combined_df.reset_index(drop=True, inplace=True)
	#df = x[1]
	df.drop(df.tail(1).index,inplace=True) 
	sv_types = df.index.tolist()
	df.insert(0, 'SV Type', sv_types)
	#df['SV Type'] = df.index
	combined_df = combined_df.append(df,ignore_index=True)

filepath = Path('4_results/full_depth/sv_intersection_agreement/supplemental/S5-uniq_prop_all_callers.csv') 
combined_df.to_csv(filepath, float_format='%.2f') 
#write_excel(EXCEL_FILEPATH,uniq_prop_all_callers,"S4")