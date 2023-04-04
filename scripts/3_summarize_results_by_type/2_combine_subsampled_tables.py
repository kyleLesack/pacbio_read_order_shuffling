import argparse
import os
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("table_10x", help="The csv file that summarizes results for data subsampled to 10X")
parser.add_argument("table_20x", help="The csv file that summarizes results for data subsampled to 20X")
parser.add_argument("table_40x", help="The csv file that summarizes results for data subsampled to 40X")
parser.add_argument("table_60x", help="The csv file that summarizes results for data subsampled to 60X")
parser.add_argument("caller", help="Variant caller")

args = parser.parse_args()

# Read csv from disk
def import_csv(filename):
	mydf = pd.read_csv(filename)
	return(mydf)

# Write results to disk
def write_results(total_df, outputdir, outputfilename,writeindex):
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	outfile = outputdir.joinpath(outputfilename)

	print("Writing to: " + str(outfile))
	total_df.to_csv(outfile, float_format="%.2f", index = writeindex)

summary_table_10x = import_csv(args.table_10x)
summary_table_20x = import_csv(args.table_20x)
summary_table_40x = import_csv(args.table_40x)
summary_table_60x = import_csv(args.table_60x)

summary_table_10x['Depth'] = "10X"
summary_table_20x['Depth'] = "20X"
summary_table_40x['Depth'] = "40X"
summary_table_60x['Depth'] = "60X"

summary_table_10x_wide  = summary_table_10x[['SV TYPE', 'Depth', 'Total calls', 'Non-intersecting proportion']]
summary_table_10x_wide["Total calls_Difference prop"] = summary_table_10x_wide["Total calls"].astype(str) + "," + summary_table_10x_wide["Non-intersecting proportion"].astype(str)
summary_table_10x_wide = summary_table_10x_wide.pivot(index='Depth', columns='SV TYPE', values=['Total calls_Difference prop'])
summary_table_20x_wide  = summary_table_20x[['SV TYPE', 'Depth', 'Total calls', 'Non-intersecting proportion']]
summary_table_20x_wide["Total calls_Difference prop"] = summary_table_20x_wide["Total calls"].astype(str) + "," + summary_table_20x_wide["Non-intersecting proportion"].astype(str)
summary_table_20x_wide = summary_table_20x_wide.pivot(index='Depth', columns='SV TYPE', values=['Total calls_Difference prop'])
summary_table_40x_wide  = summary_table_40x[['SV TYPE', 'Depth', 'Total calls', 'Non-intersecting proportion']]
summary_table_40x_wide["Total calls_Difference prop"] = summary_table_40x_wide["Total calls"].astype(str) + "," + summary_table_40x_wide["Non-intersecting proportion"].astype(str)
summary_table_40x_wide = summary_table_40x_wide.pivot(index='Depth', columns='SV TYPE', values=['Total calls_Difference prop'])
summary_table_60x_wide  = summary_table_60x[['SV TYPE', 'Depth', 'Total calls', 'Non-intersecting proportion']]
summary_table_60x_wide["Total calls_Difference prop"] = summary_table_60x_wide["Total calls"].astype(str) + "," + summary_table_60x_wide["Non-intersecting proportion"].astype(str)
summary_table_60x_wide = summary_table_60x_wide.pivot(index='Depth', columns='SV TYPE', values=['Total calls_Difference prop'])

total_df = pd.concat([summary_table_10x,summary_table_20x,summary_table_40x, summary_table_60x])
total_df_wide = pd.concat([summary_table_10x_wide,summary_table_20x_wide,summary_table_40x_wide, summary_table_60x_wide])
total_df_wide.columns = total_df_wide.columns.droplevel(0)
total_df_wide.reset_index(inplace=True)

total_df = total_df[['SV TYPE', 'Depth', 'Intersection', 'Non-intersecting', 'Non-intersecting proportion']]
outputpath = Path(args.table_10x.replace("10X", "combined"))
outputfilename = outputpath.name
outputdir = outputpath.parent.absolute().joinpath("combined_long_table") # Write long format table to disk
write_results(total_df, outputdir, outputfilename, False)
outputdir = outputpath.parent.absolute().joinpath("combined_wide_table") # Write wide format table to disk
write_results(total_df_wide, outputdir, outputfilename, False)
