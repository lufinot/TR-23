import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from dbscan1d.core import DBSCAN1D

import argparse
import os
import logging
import warnings

def calculate_wilcoxon_pvals(df: pd.DataFrame) -> np.ndarray:
    """
    Calculate Wilcoxon p-values for each column in a dataframe.
    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        pvals: Array of p-values for each region.
    """
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("once")
        pvals = wilcoxon(df, nan_policy='omit', zero_method='pratt', axis=0)[1]
        for warning in w:
            if str(warning.message) == 'Exact p-value calculation does not work if there are zeros. Switching to normal approximation.':
                logging.warning('Exact p-value calculation switched to normal approximation due to presence of zeros.')
            elif str(warning.message) == 'Sample size too small for normal approximation.':
                logging.warning('Sample size is too small for normal approximation.')
    _, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    return pvals_corrected



def process_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove columns with all zeros and sample_id column if present.
    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        df: Processed dataframe.
    """
    if 'sample_id' in df.columns:
        df.drop(columns=['sample_id'], axis=1, inplace=True)
    df = df.loc[:, (df.fillna(0) != 0).any(axis=0)]
    return df


# Need to add motif extraction to GenotypeGetter & Remove this
def add_motif_info(df: pd.DataFrame) -> pd.DataFrame:
    """ 
    Add motif information to the dataframe, from the locus_structures.csv file.
    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        df: Dataframe with motif information added.
    """
    loc_to_motif = pd.read_csv('data/other/locus_structures.csv')
    loc_to_motif['ReferenceRegion'] = loc_to_motif['ReferenceRegion'].str.lstrip('chr')
    df = df.merge(loc_to_motif[['ReferenceRegion', 'LocusStructure']], on='ReferenceRegion', how='left')
    return df


def cluster_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Call cluster function to get cluster features for each locus. 
    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        cluster_df: With Cols ['num_clusters', 'cluster_means', 'cluster_sds', 'out4', 'out6']
    """
    result = df.apply(lambda x: cluster_and_outliers(x), axis=0).T
    result.columns = ['num_clusters', 'cluster_means', 'cluster_sds', 'out4', 'out6']
    return result



   


def cluster_and_outliers(x: pd.Series) -> list:
    x = x.dropna()
    sx = x.div(x.std())
    sx = sx.values
    db = DBSCAN1D(eps=0.3, min_samples=20)
    labels = db.fit_predict(sx)

    outliers = sx[labels == -1]

    # get counts of outliers above 4 and 6 sd
    out4 = len(outliers[outliers > 4])
    out6 = len(outliers[outliers > 6])

    # get mean and sd of clusters
    means = []
    sds = []
    for label in np.unique(labels)[1:]:
        means.append(np.mean(x[labels == label]))
        sds.append(np.std(x[labels == label]))

    # remove from means and sd list the cluster with the mean closest to 0
    df = pd.DataFrame({
        'means': means,
        'sds': sds
    })
    df['abs_means'] = df['means'].abs()
    df = df.sort_values('abs_means')
    df = df.iloc[1:, :]
 
    means = df['means'].values.tolist()  # Convert np.array to list
    sds = df['sds'].values.tolist()  # Convert np.array to list

    return [len(means), means, sds, out4, out6]  # Return a list instead of a tuple

    


def process_and_extract_features(input_df_path: str) -> pd.DataFrame:
    """
    Function to call other functions to process the dataframe and extract features.
    Args:
        input_df_path: Path to the dataframe with rows as samples and cols as regions.
    Returns:
        features_df: Dataframe with loci as rows, with features as cols.
    """
    # Read in difference dataframe
    df = process_df(df)

    features_df = args.feats or pd.DataFrame({'ReferenceRegion': df.columns})

    features_df = add_motif_info(features_df)

    features_df['wilcox_pvals'] = calculate_wilcoxon_pvals(df)

    features_df = features_df.concat([features_df, cluster_features(df)], axis=1)


    # get normalized non-zero proportion
    features_df['prop_nonzero'] = df.astype(bool).sum(axis=0)/(df.shape[0]+(2*np.sqrt(df.shape[0])))

    # get standard deviation
    features_df['std'] = df.std(axis=0)

    return features_df


def init_argparse():
    parser = argparse.ArgumentParser(description='Calculate the p-values from tidied data.')
    parser.add_argument('--tidat', required=True, help='Location of the tidied data')
    parser.add_argument('--outdir', default='', help='Output directory for the pvals. (default: script running directory)')
    parser.add_argument('--name', help='Prefix for output file (default same as tidied dat)')
    parser.add_argument('--feats', help='Location of the features file (default: creates new features)')
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
