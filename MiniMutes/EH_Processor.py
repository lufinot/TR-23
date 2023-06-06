### Script to process the EH data which has been converted to an ndjson file. 
### Input: folder containing the ndjson files
###        a csv file containing the list of files to be processed

import pandas as pd
import sys

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


def main(ExpansionHunterData, manifest):
    # Load and diff the data
    mani = pd.read_csv(manifest)
    mani = mani.head(3)
    # Create an empty DataFrame 
    dfs = []
    for index, row in mani.iterrows():
        dfs.append(load_and_diff(ExpansionHunterData, row))
        print("Processed {} out of {}".format(index+1, len(mani)))
    df = pd.concat(dfs, axis=0)
    df.to_csv('EHPoutput.csv')

    # Loop through the manifest and load and diff the data
    print("Saved to EHPoutput.csv")



if __name__ == "__main__":
   mani, loc = sys.argv[1:3]
   print(main(loc, mani))


