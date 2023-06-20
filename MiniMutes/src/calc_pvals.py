import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import diptest
import argparse
import os
import logging
import warnings


def calculate_wilcoxon_pvals(df: pd.DataFrame) -> np.ndarray:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("once")
        pvals = wilcoxon(df, nan_policy='omit', axis=0)[1]
        for warning in w:
            if str(warning.message) == 'Exact p-value calculation does not work if there are zeros. Switching to normal approximation.':
                logging.warning('Exact p-value calculation switched to normal approximation due to presence of zeros.')
            elif str(warning.message) == 'Sample size too small for normal approximation.':
                logging.warning('Sample size is too small for normal approximation.')
    _, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    return pvals_corrected


def calculate_dip_pvals(df: pd.DataFrame) -> np.ndarray:
    pvals = df.apply(lambda x: diptest.diptest(x.values)[1], axis=0)
    _, pvals_corrected, _, _ = multipletests(pvals.values, alpha=0.05, method='fdr_bh')
    return pvals_corrected


def process_df(df: pd.DataFrame) -> pd.DataFrame:
    if 'sample_id' in df.columns:
        df.drop(columns=['sample_id'], axis=1, inplace=True)
    df = df.loc[:, (df.fillna(0) != 0).any(axis=0)]
    return df


def add_motif_info(df: pd.DataFrame) -> pd.DataFrame:
    loc_to_motif = pd.read_csv('data/other/locus_structures.csv')
    loc_to_motif['ReferenceRegion'] = loc_to_motif['ReferenceRegion'].str.lstrip('chr')
    df = df.merge(loc_to_motif[['ReferenceRegion', 'LocusStructure']], on='ReferenceRegion', how='left')
    return df


def process_and_get_pvals(input_df_path: str) -> pd.DataFrame:
    
    df = process_df(df)

    pvals_df = pd.DataFrame({'ReferenceRegion': df.columns})
    pvals_df['wilcox_pvals'] = calculate_wilcoxon_pvals(df)
    pvals_df['dip_pvals'] = calculate_dip_pvals(df)

    pvals_df = add_motif_info(pvals_df)

    # get normalized non-zero proportion
    pvals_df['prop_nonzero'] = df.astype(bool).sum(axis=0)/(df.shape[0]+(2*np.sqrt(df.shape[0])))

    # get standard deviation
    pvals_df['std'] = df.std(axis=0)

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
    pvals_df = process_and_get_pvals(df)
    output_name = args.name if args.name else os.path.basename(args.tidat).split('.')[0]
    output_path = os.path.join(args.outdir, f"{output_name}_pvals.csv")

    pvals_df.to_csv(output_path, index=False)

    print(f"P-values saved to {output_path}")


if __name__ == "__main__":
    main()
