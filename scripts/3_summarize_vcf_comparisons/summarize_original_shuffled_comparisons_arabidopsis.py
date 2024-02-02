import argparse
import pandas

ALL_STRAINS = ["1254", "6021", "6024", "9470"]

parser = argparse.ArgumentParser()
parser.add_argument("input_csv", help=r"Input csv file format. Placeholder text should be used for the strain (STRAIN) and replicate (REPLICATE)") # 4_variant_call_comparisons/subsampled/20X/pbmm2/pbsv/STRAIN/original_shuffled_comparisons/shuffledREPLICATE/results_proportions.csv
parser.add_argument("output_dir", help="Output directory")
parser.add_argument("--num_shuffled", default=5, type=int, help="Number of shuffled files")

args = parser.parse_args()

def parse_csv_file(csv_file, strain):
	with open(csv_file) as f:
		for line in f:
			pass
	last_line = line.rstrip().split(",")

	if "Total" not in last_line: # Last line should contain the totals. Use this as a sanity check.
		print("WARNING: Totals not in last line")
	last_line[0] = strain

	return last_line

csv_totals_list = []

for strain in ALL_STRAINS:
	for replicate in range(1, args.num_shuffled + 1):
		input_file = args.input_csv
		input_file = input_file.replace("STRAIN", strain).replace("REPLICATE", str(replicate))
		csv_totals_list.append(parse_csv_file(input_file, strain))

# Create data frame with Jaccard statistics and write to disk
csv_totals_df = pandas.DataFrame(csv_totals_list, columns=['Strain','JaccardIndex','JaccardDistance','Symmetric Difference']) # Create data frame with distances and differences
csv_totals_df.drop('Symmetric Difference', axis=1, inplace=True)
csv_totals_file = args.output_dir + "/csv_totals_all_strains_replicates_jaccard.csv"
csv_totals_df.to_csv(csv_totals_file, index = False, float_format='%f')

# Create data frame with symmetric differences and write to disk
csv_totals_df = pandas.DataFrame(csv_totals_list, columns=['Strain','JaccardIndex','JaccardDistance','Symmetric Difference']) # Create data frame with distances and differences
csv_totals_df.drop(['JaccardIndex', 'JaccardDistance'], axis=1, inplace=True)
csv_totals_file = args.output_dir + "/csv_totals_all_strains_replicates_differences.csv"
csv_totals_df.to_csv(csv_totals_file, index = False, float_format='%f')
