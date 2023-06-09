# calc_pvals.py
import pandas as pd
from scipy.stats import wilcoxon
import argparse
import os

def get_pvals(df) -> pd.DataFrame:
    # Remove the sample_id column if it exists
    if 'sample_id' in df.columns:
        df.drop(columns=['sample_id'], axis=1, inplace=True)

    # Deal with NaNs and 0 cols
    dfVals = df.fillna(0)
    dfVals = df.loc[:, (dfVals != 0).any(axis=0)]

    # Get the p-values
    pvals = wilcoxon(dfVals.iloc[:, 1:], nan_policy='omit', axis = 0)[1]
    pvals = pd.DataFrame({'pval': pvals}, index=dfVals.columns[1:])
    return pvals

def init_argparse():
    parser = argparse.ArgumentParser(description='Calculate the p-values from tidied data.')
    parser.add_argument('--tidat', required=True, help='Location of the tidied data')
    return parser

def main(args=None):
    parser = init_argparse()
    args = parser.parse_args()

    if not os.path.exists(args.tidat):
        print(f"Error: {args.tidat} does not exist.")
        return
    df = pd.read_csv(args.tidat, index_col=0)

    # Process the tidied data and save pvals to a csv
    pvals = get_pvals(df)
    pvals.to_csv("{prefix}_pvals.csv".format(prefix=os.path.splitext(args.tidat)[0]))
    print("P-values saved to {prefix}_pvals.csv".format(prefix=os.path.splitext(args.tidat)[0]))


if __name__ == "__main__":
    main()
