import pandas as pd
import argparse
import logging
import os

def load_dataframe(path):
    try:
        df = pd.read_csv(path, index_col=0)
        if df.empty:
            logging.error(f"{path} is empty. ")
            return None
        return df
    except (FileNotFoundError, ValueError) as err:
        logging.error(f"{err}")
        return None

def process_dataframe(df_case, df_control):
    # Subtract genotype columns
    if df_case.shape != df_control.shape:
        logging.error(f"Case and control dataframes have different shapes. Can't perform subtraction.")
        return None

    df_diff = df_case.subtract(df_control)

    return df_diff

def init_argparse():
    parser = argparse.ArgumentParser(
        usage="%(prog)s --case <case csv> --control <control csv> --outdir <output directory>",
        description='Differences the genotype columns of the case and control dataframes.'
    )
    parser.add_argument('--case', required=True, help='CSV file containing case genotypes')
    parser.add_argument('--control', required=True, help='CSV file containing control genotypes')
    parser.add_argument('--name', required=True, help='Output path for the tidied diff. (default: script running directory)')
    return parser

def main():
    parser = init_argparse()
    args = parser.parse_args()

    # Load the dataframes
    df_case = load_dataframe(args.case)
    df_control = load_dataframe(args.control)

    # If any dataframe could not be loaded, exit
    if df_case is None or df_control is None:
        exit()

    # Process the dataframes ignoring first column
    df_diff = process_dataframe(df_case.iloc[:,1:], df_control.iloc[:,1:])

    df_diff.insert(0, 'sample_id', df_case['sample_id'])
    # If dataframe couldn't be processed, exit
    if df_diff is None:
        exit()


    df_diff.to_csv(args.name + '.csv')
    print(f"Differenced df saved to {args.name}.csv")

if __name__ == "__main__":
    main()
