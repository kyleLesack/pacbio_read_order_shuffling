import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-p','--pbsv_files', nargs='+', help='<Required> List of pbsv csv files to summarize', required=True)
parser.add_argument('-s','--sniffles_files', nargs='+', help='<Required> List of Sniffles csv files to summarize', required=True)
parser.add_argument('-v','--svim_files', nargs='+', help='<Required> List of SVIM csv files to summarize', required=True)
parser.add_argument('-o','--output_dir', help='<Required> Output directory', required=True)
parser.add_argument('-n','--strain', help='<Required> Strain name', required=True)
parser.add_argument('--no_meta', action=argparse.BooleanOptionalAction)

args = parser.parse_args()
pbsv_csv_files = args.pbsv_files
sniffles_csv_files = args.sniffles_files
svim_csv_files = args.svim_files

def pop_std(x): # Function to calculate population standard deviation
	return x.std(ddof=0)

# Read data from csv file
def parse_csv_files(csv_files,caller):
	csv_lines_list = []
	for file in csv_files:
		with open(file) as f:
			sv_stats = f.readlines()
			csv_header = sv_stats[0].rstrip().split(",")
			for line in sv_stats[1:]:
				line_split = line.rstrip().split(",")
				csv_lines_list.append(line_split)

	sv_agreement_stats_df = pd.DataFrame(csv_lines_list, columns = csv_header)
	sv_agreement_stats_df["Jaccard Index"] = pd.to_numeric(sv_agreement_stats_df["Jaccard Index"])
	sv_agreement_stats_df["Jaccard Distance"] = pd.to_numeric(sv_agreement_stats_df["Jaccard Distance"])
	sv_agreement_stats_df["Symmetric Difference"] = pd.to_numeric(sv_agreement_stats_df["Symmetric Difference"])

	return(sv_agreement_stats_df)

# Get Jaccard distances in wide format
def get_jaccard_dist_wide(caller_agreement_stats_df, caller):
	caller_agreement_stats_jaccard_distance_df = caller_agreement_stats_df[['SV Type', 'Jaccard Distance']]
	caller_agreement_stats_jaccard_distance_df = caller_agreement_stats_jaccard_distance_df.groupby('SV Type').agg(["mean", pop_std]) # Creates new df where the rows are the SV Types and columns are the mean and standard deviation of the values
	caller_agreement_stats_jaccard_distance_df = caller_agreement_stats_jaccard_distance_df.loc[:, 'Jaccard Distance']
	caller_agreement_stats_jaccard_distance_df['mean'] = caller_agreement_stats_jaccard_distance_df['mean'].apply(lambda x: '{0:.3f}'.format(x)) # Convert to string to allow mean + standard dev to be stored in same column
	caller_agreement_stats_jaccard_distance_df['pop_std'] = caller_agreement_stats_jaccard_distance_df['pop_std'].apply(lambda x: '{0:.3f}'.format(x)) # Convert to string to allow mean + standard dev to be stored in same column
	caller_agreement_stats_jaccard_distance_df["Jaccard Distance"] = caller_agreement_stats_jaccard_distance_df['mean'] + " ± " + caller_agreement_stats_jaccard_distance_df['pop_std'] # Create new col with mean +- std
	caller_agreement_stats_jaccard_distance_df['Caller'] = caller
	caller_agreement_stats_jaccard_distance_df = caller_agreement_stats_jaccard_distance_df.reset_index()
	jaccard_dist_caller_df = caller_agreement_stats_jaccard_distance_df[["SV Type", "Caller", "Jaccard Distance"]]
	jaccard_dist_caller_df = jaccard_dist_caller_df.pivot(index='Caller', columns='SV Type', values='Jaccard Distance')
	return(jaccard_dist_caller_df)

# Import data for each caller and convert to wide format
pbsv_agreement_stats_df = parse_csv_files(pbsv_csv_files, "pbsv")
jaccard_dist_pbsv_df = get_jaccard_dist_wide(pbsv_agreement_stats_df, "pbsv")
sniffles_agreement_stats_df = parse_csv_files(sniffles_csv_files, "sniffles")
jaccard_dist_sniffles_df = get_jaccard_dist_wide(sniffles_agreement_stats_df, "sniffles")
svim_agreement_stats_df = parse_csv_files(svim_csv_files, "svim")
jaccard_dist_svim_df = get_jaccard_dist_wide(svim_agreement_stats_df, "svim")
combined_jaccard_dist_df = pd.concat([jaccard_dist_pbsv_df, jaccard_dist_sniffles_df, jaccard_dist_svim_df])
combined_jaccard_dist_df["Strain"] = args.strain
Path(args.output_dir).mkdir(parents=True, exist_ok=True)

if args.no_meta:
	combined_jaccard_dist_file = args.output_dir + "/" + args.strain + "_results_no_meta.csv"
	combined_jaccard_dist_df.to_csv(combined_jaccard_dist_file, index=True)
else:
	combined_jaccard_dist_file = args.output_dir + "/" + args.strain + "_results.csv"
	combined_jaccard_dist_df.to_csv(combined_jaccard_dist_file, index=True)
