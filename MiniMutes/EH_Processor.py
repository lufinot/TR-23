### Script to process the EH data which has been converted to an ndjson file. 
### Input: folder containing the ndjson files
###        a csv file containing the list of files to be processed

import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import sys

# Split the genotype column into two numeric columns
def process_dataframe(df):
    # Apply add_slash to the 'genotype' column
    df['genotype'] = df['genotype'].astype('string').apply(add_slash)

    # Split the 'genotype' column into two
    df[['genotype1', 'genotype2']] = df['genotype'].str.split('/', expand=True)

    # Convert genotypes to numeric, errors='coerce' will convert invalid parsing to NaN
    df[['genotype1', 'genotype2']] = df[['genotype1', 'genotype2']].apply(pd.to_numeric, errors='coerce')

    return df


# Load the data, check for errors, and return the DataFrame
def load_dataframe(path):
    try:
        df = pd.read_json(path, lines=True)
        # Check if 'genotype' column exists
        if 'genotype' not in df.columns:
            print(f"Error: {path} does not contain a 'genotype' column.", file=sys.stderr)
            return None
        # Check if DataFrame is empty
        if df.empty:
            print(f"Error: {path} is empty.", file=sys.stderr)
            return None
        return df
    except (FileNotFoundError, IsADirectoryError, ValueError) as err:
        print(f"Error: {path}: {err.strerror}", file=sys.stderr)
        return None


# if there is no slash in the string, make it 'orginal/original'
def add_slash(string):
    if pd.isna(string):
        return string
    elif '/' not in string:
        return '{}/{}'.format(string, string)
    else:
        return string
    
    
# Load the data and subtract the genotype columns
# loc: folder with ndjson files
# sample: row from the mainfest
def load_and_diff(sample, loc):
    print(type(sample))  # Add this line
    print(sample)        # And this line

    path_format = loc + '/{id}.ndjson'
    print('Processing {}'.format(sample['icgc_donor_id']))

    # Load the dataframes
    dfn_raw = load_dataframe(path_format.format(id=sample['control_object_id']))
    dft_raw = load_dataframe(path_format.format(id=sample['case_object_id']))

    # If any dataframe could not be loaded, return None
    if dfn_raw is None or dft_raw is None:
        return None

    # Process the dataframes
    dfn = process_dataframe(dfn_raw)
    dft = process_dataframe(dft_raw)

    # Subtract genotype1 and genotype2 columns in both dataframes
    diff1 = dfn['genotype1'].subtract(dft['genotype1'])
    diff2 = dfn['genotype2'].subtract(dft['genotype2'])

    # Create a new DataFrame from these Series
    row_names = ['0' + sample['icgc_donor_id'], '1' + sample['icgc_donor_id']]
    diff_df = pd.DataFrame([diff1, diff2], index = row_names)
    diff_df.columns = dfn['region']
    diff_df.insert(0, 'sample_id', sample['icgc_donor_id'])

    return diff_df


def get_pvals(df) -> pd.DataFrame:
    # Remove the sample_id column if it exists
    if 'sample_id' in df.columns:
        df = df.drop(columns=['sample_id'], axis=1, inplace=True)
    
    # Deal with NaNs and 0 cols
    dfVals = df.fillna(0)
    dfVals = df.loc[:, (dfVals != 0).any(axis=0)]

    # Get the p-values
    pvals = wilcoxon(dfVals.iloc[:, 1:], nan_policy='omit', axis = 0)[1]
    pvals = pd.DataFrame({'pval': pvals}, index=dfVals.columns[1:])
    return pvals

def init_argparse():
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description='Process the EH data.')
    parser.add_argument('--EHD', required=True, help='Directory w/ ndjson files')
    parser.add_argument('--manifest', required=True, help='CSV of paired files to be processed, with columns: icgc_donor_id, control_object_id, case_object_id, sex')
    parser.add_argument('--dis', help='Name of the disease (default: Name of --EHD arg)')
    parser.add_argument('--step', default=0, type=int, choices=[0, 1], help='Step to process the data from: {0: raw ndjson, 1: tidied data} (default 0)')
    parser.add_argument('--tidat', help='Location of the tidied data (If step is 1)')
    
    return parser

def main():
    parser = init_argparse()
    args = parser.parse_args()

    # Set the disease name if it was not given
    if args.dis is None:
        args.dis = os.path.basename(args.EHD)

    # Depending on the step, do the appropriate processing
    if args.step == 0:
        # Load and diff the data
        mani = pd.read_csv(args.manifest)

        dfs = mani.apply(lambda row: load_and_diff(row, args.EHD), axis=1).tolist()
        dfs = [df for df in dfs if df is not None]
        df = pd.concat(dfs, axis=0)
        df.to_csv("{disease}_res.csv".format(disease=args.dis))

        print("Tidied df saved to {disease}.csv".format(disease=args.dis))
    elif args.step == 1:
        # Check if the tidied_data argument was given
        if args.tidat is None:
            print("Error: the --tidat argument is necessary when --step is 1.")
            return
        if not os.path.exists(args.tidat):
            print(f"Error: {args.tidat} does not exist.")
            return
        df = pd.read_csv(args.tidat, index_col=0)
    
    # Process the tidied data and save pvals to a csv
    pvals = get_pvals(df)
    pvals.to_csv("{disease}_pvals.csv".format(disease=args.dis))
    print("P-values saved to {disease}_pvals.csv".format(disease=args.dis))


if __name__ == "__main__":
    main()