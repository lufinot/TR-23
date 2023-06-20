### Script to process the EH data which has been converted to an ndjson file. 
### Input: folder containing the ndjson files
###        a csv file containing the list of files to be processed

import os
import sys
import argparse
import pandas as pd
from calc_pvals import get_pvals
import logging 

LOG_LEVEL = os.getenv('LOG_LEVEL') or 'info'

log_dict = {'debug':logging.DEBUG, 'info':logging.INFO , 'warning':logging.WARNING, 
            'error':logging.ERROR, 'critical':logging.CRITICAL}

log_level = log_dict.get(LOG_LEVEL.lower(),logging.INFO)

logging.basicConfig(filename='filename.log', level=log_level,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

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
            logging.error(f"{path} does not contain a 'genotype' column. ")
            return None
        # Check if DataFrame is empty
        if df.empty:
            logging.error(f"{path} is empty. ")
            return None
        return df
    except (FileNotFoundError, IsADirectoryError, ValueError) as err:
        logging.error(f"{err}")
        return None


# if there is no slash in the string, make it 'orginal/original'
def add_slash(string):
    if pd.isna(string):
        return string
    elif '/' not in string:
        return f'{string}/NaN'
    else:
        return string
    
    
# Load the case data and format the DataFrame
# loc: folder with ndjson files
# sample: row from the mainfest
def load_and_case(sample, loc):
    path_format = loc +'/{id}.ndjson'
    logging.info(f"Loading {sample['icgc_donor_id']}...")

    # Load the dataframe
    dft_raw = load_dataframe(path_format.format(id=sample['case_object_id']))

    # If dataframe could not be loaded, return None
    if dft_raw is None:
        logging.warning(f"{sample['icgc_donor_id']} could not be loaded.")
        return None

    # Process the dataframe
    dft = process_dataframe(dft_raw)

    totNaNs = dft.isna().any(axis=1).sum()
    # log if more than 30% of the rows have NaN values
    if totNaNs > 0.3 * len(dft):
        logging.warning(f"{sample['icgc_donor_id']} has more than 30% missing values.")

    # Create a new DataFrame from these Series
    row_names = ['0' + sample['icgc_donor_id'], '1' + sample['icgc_donor_id']]
    case_df = pd.DataFrame([dft['genotype1'], dft['genotype2']], index = row_names)
    case_df.columns = dft['region']
    case_df.insert(0, 'sample_id', sample['icgc_donor_id'])

    return case_df

def init_argparse():
    parser = argparse.ArgumentParser(
        usage="%(prog)s --EHD <EHD directory> --manifest <manifest file> [other options]]",
        description='Process Expansion Hunter output reformatted to ndjson files into a tidy DataFrame. Also calculates p-values.')
    parser.add_argument('--EHD', required=True, help='Directory w/ ndjson files')
    parser.add_argument('--manifest', required=True, help='CSV of paired files to be processed, with columns: icgc_donor_id, control_object_id, case_object_id, sex')
    parser.add_argument('--disease', help='Name of the disease (default: Name of --EHD arg)')
    parser.add_argument('--pvals', default=True, type=bool, help='Calculate p-values? (default True)')
    parser.add_argument('--outdir', default='', help='Output directory for the tidied data and pvals. (default: script running directory)')
    return parser

def main():
    parser = init_argparse()
    args = parser.parse_args()

    # Set the disease name if it was not given
    args.disease = args.disease or os.path.basename(args.EHD)

    # Load and diff the data
    mani = pd.read_csv(args.manifest)
    dfs = mani.apply(lambda row: load_and_case(row, args.EHD), axis=1).tolist()

    dfs = [df for df in dfs if df is not None]
    df = pd.concat(dfs, axis=0)
    
    output_dir = args.outdir if args.outdir else "."
    tidied_file = os.path.join(output_dir, f"{args.disease}_tidy.csv")
    df.to_csv(tidied_file)
    print(f"Tidied df saved to {tidied_file}")


if __name__ == "__main__":
    main()
