# calc_pvals.py
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import argparse
import os
import logging
import warnings

def get_pvals(df) -> pd.DataFrame:
    # Remove the sample_id column if it exists
    if 'sample_id' in df.columns:
        df.drop(columns=['sample_id'], axis=1, inplace=True)

    # Deal with NaNs and 0 cols
    dfVals = df.fillna(0)
    dfVals = df.loc[:, (dfVals != 0).any(axis=0)]

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("once") 
        pvals = wilcoxon(dfVals.iloc[:, 1:], nan_policy='omit', axis = 0)[1]
        for warning in w:
            if str(warning.message) == 'Exact p-value calculation does not work if there are zeros. Switching to normal approximation.':
                logging.warning('Exact p-value calculation switched to normal approximation due to presence of zeros.')
            elif str(warning.message) == 'Sample size too small for normal approximation.':
                logging.warning('Sample size is too small for normal approximation.')
    
    # Apply FDR correction
    rejected, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    
    pvals_df = pd.DataFrame({'pval': pvals, 'pval_corrected': pvals_corrected, 'rejected': rejected}, index=dfVals.columns[1:])
    
    
    return pvals_df


def init_argparse():
    parser = argparse.ArgumentParser(description='Calculate the p-values from tidied data.')
    parser.add_argument('--tidat', required=True, help='Location of the tidied data')
    parser.add_argument('--outdir', default='', help='Output directory for the pvals. (default: script running directory)')
    parser.add_argument('--name', help='Prefix for output file (default same as tidied dat)')
    return parser

def main(args=None):
    logging.basicConfig(filename='calc_pvals.log', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    parser = init_argparse()
    args = parser.parse_args()

    if not os.path.exists(args.tidat):
        print(f"Error: {args.tidat} does not exist.")
        return
    df = pd.read_csv(args.tidat, index_col=0)
    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.name or os.path.splitext(args.tidat)[0]
    path = f'{args.outdir}/{prefix}.csv'
    # Process the tidied data and save pvals to a csv
    pvals = get_pvals(df)
    pvals.to_csv(path)
    print(f'P-values saved to {path}')


if __name__ == "__main__":
    main()
