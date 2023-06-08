### Script to process the EH data which has been converted to an ndjson file. 
### Input: folder containing the ndjson files
###        a csv file containing the list of files to be processed

import pandas as pd
import os
import argparse

# Load the data and split the genotype column into two numeric columns
def process_dataframe(path):
    # Load the data
    df = pd.read_json(path, lines=True)

    # Apply add_slash to the 'genotype' column
    df['genotype'] = df['genotype'].astype('string').apply(add_slash)

    # Split the 'genotype' column into two
    df[['genotype1', 'genotype2']] = df['genotype'].str.split('/', expand=True)

    # Convert genotypes to numeric, errors='coerce' will convert invalid parsing to NaN
    df[['genotype1', 'genotype2']] = df[['genotype1', 'genotype2']].apply(pd.to_numeric, errors='coerce')

    return df



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
def load_and_diff(loc, sample):
    path_format = loc + '/{id}.ndjson'
    print('Processing {}'.format(sample['icgc_donor_id']))
    # Use the helper function to process the dataframes
    dfn = process_dataframe(path_format.format(id=sample['control_object_id']))
    dft = process_dataframe(path_format.format(id=sample['case_object_id']))

    # Subtract genotype1 and genotype2 columns in both dataframes
    diff1 = dfn['genotype1'].subtract(dft['genotype1'])
    diff2 = dfn['genotype2'].subtract(dft['genotype2'])

    # Create a new DataFrame from these Series
    row_names = ['0' + sample['icgc_donor_id'], '1' + sample['icgc_donor_id']]
    diff_df = pd.DataFrame([diff1, diff2], index = row_names)
    diff_df.columns = dfn['region']
    diff_df.insert(0, 'sample_id', sample['icgc_donor_id'])

    return diff_df


def init_argparse():
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description='Process the EH data.')
    parser.add_argument('--ExpansionHunterData', required=True, help='Directory w/ ndjson files')
    parser.add_argument('--manifest', required=True, help='Paired files to be processed')
    parser.add_argument('--disease_name', required=True, help='Name of the disease')
    parser.add_argument('--step', default=0, type=int, choices=[0, 1], help='Step to process the data from: {0: raw ndjson, 1: tidied data} (default 0)')
    parser.add_argument('--tidied_data', help='Location of the tidied data (If step is 1)')
    
    return parser

def main():
    parser = init_argparse()
    args = parser.parse_args()

    # Check if the files exist
    if not os.path.exists(args.ExpansionHunterData):
        print(f"Error: {args.ExpansionHunterData} does not exist.")
        return
    if not os.path.exists(args.manifest):
        print(f"Error: {args.manifest} does not exist.")
        return

    # Depending on the step, do the appropriate processing
    if args.step == 0:
        # Load and diff the data
        mani = pd.read_csv(args.manifest)
        # Create an empty DataFrame 
        dfs = []
        for index, row in mani.iterrows():
            dfs.append(load_and_diff(args.ExpansionHunterData, row))
        df = pd.concat(dfs, axis=0)
        df.to_csv("{disease}_res.csv".format(disease=args.disease_name))

        # Loop through the manifest and load and diff the data
        print("Saved to {disease}.csv".format(disease=args.disease_name))
    elif args.step == 1:
        # Check if the tidied_data argument was given
        if args.tidied_data is None:
            print("Error: the --tidied_data argument is necessary when --step is 1.")
            return
        if not os.path.exists(args.tidied_data):
            print(f"Error: {args.tidied_data} does not exist.")
            return

        # Here you would include the code for processing tidied data
        pass

if __name__ == "__main__":
    main()