import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_files', nargs='+', help='<Required> List of input files generated from the original and shuffled per strain agreement results', required=True)
parser.add_argument('-o','--output_file', help='<Required> Output file', required=True)

args = parser.parse_args()

# Parse CSV file
def parse_csv_files(csv_files):
	csv_lines_list = []
	for file in csv_files:
		with open(file) as f:
			csv_lines = f.readlines()
			csv_header = csv_lines[0].rstrip().split(",")

			for line in csv_lines[1:]:
				line_split = line.rstrip().split(",")
				csv_lines_list.append(line_split)

	csv_lines_df = pd.DataFrame(csv_lines_list, columns = csv_header)
	cols = list(csv_lines_df.columns.values) # Get list of column names in order to sort columns alphabetically, so that SVIM duplication columns aren't last
	cols.sort()
	csv_lines_df = csv_lines_df.reindex(columns=cols)
	cols_to_move = ['Strain', 'Caller']
	csv_lines_df = csv_lines_df[ cols_to_move + [ col for col in csv_lines_df.columns if col not in cols_to_move ] ] # https://stackoverflow.com/a/56479671

	return(csv_lines_df)

csv_lines_df = parse_csv_files(args.input_files)

output_file = args.output_file
path = Path(output_file)
output_dir = path.parent.absolute()
Path(output_dir).mkdir(parents=True, exist_ok=True)
csv_lines_df.to_csv(path, index=False)
