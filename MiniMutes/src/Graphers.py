import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


### GENOTYPE GRAPHING FUNCTIONS ###
def graphGenotypes(loc, cancer):
    case = pd.read_csv(f'data/genotypes/{cancer}_case_tidy.csv', usecols=[loc])
    control = pd.read_csv(f'data/genotypes/{cancer}_control_tidy.csv', usecols=[loc])

    #graph the case and control genotypes on the same hist
    plt.hist(case[loc], bins=100, alpha=0.5, label='case', color='red')
    plt.hist(control[loc], bins=100, alpha=0.5, label='control', color='blue')
    plt.title(f'{cancer} ; {loc} Distribution')
    plt.xlabel('Genotype')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    
    plt.show()


def graphGenotypesLociHelper(loc, cancer, ax):
    case = pd.read_csv(f'data/genotypes/{cancer}_case_tidy.csv', usecols=[loc])
    control = pd.read_csv(f'data/genotypes/{cancer}_control_tidy.csv', usecols=[loc])

    #graph the case and control genotypes on the same hist
    ax.hist(case[loc], bins=100, alpha=0.5, label='case', color='red')
    ax.hist(control[loc], bins=100, alpha=0.5, label='control', color='blue')
    ax.set_title(f'{loc} Distribution')
    ax.set_xlabel('Genotype')
    ax.set_ylabel('Frequency')
    ax.legend(loc='upper right')

def load_graphGenotypesLoci(loci, cancer):
    num_loci = len(loci)
    num_cols = 4
    num_rows = math.ceil(num_loci / num_cols)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(30, num_rows*5)) # Adjust the figure size as per your requirement
    axs = axs.ravel() # flatten the array of axes

    for idx, loc in enumerate(loci):
        ax = axs[idx]
        graphGenotypesLociHelper(loc, cancer, ax)

    plt.title(f'{cancer} Distribution for {num_loci} Loci')
    plt.tight_layout()  # To prevent overlap of subplots
    plt.show()

def calculate_data(loci, case_df, control_df):
    data_dict = {}
    max_diff = 0
    
    for i, loc in enumerate(loci):
        # Extract data for current loci and drop null rows
        df = pd.concat([control_df[[loc]], case_df[[loc]]], axis=1)
        df.columns = ['control', 'case']
        df.dropna(inplace=True)

        # Melt the dataframe to long format for plotting
        df_melt = df.melt(value_name='Repeat Length', var_name='group')
        bins = np.arange(df_melt['Repeat Length'].min() - 0.5, df_melt['Repeat Length'].max() + 0.5)

        
        # Get histogram data for case and control groups
        hist_case, case_bins = np.histogram(df['case'], bins=bins)
        hist_control, _ = np.histogram(df['control'], bins=bins)

        # Calculate the difference in bin counts
        diff = hist_case - hist_control

        # create a convolution of the difference array with a gaussian kernel
        # to smooth out the line plot
        conv_diff = np.convolve(diff, np.ones(9)/9, mode='same')

        # if magnitude of max of conv_diff is greater than max_diff, update max_diff
        if np.abs(conv_diff).max() > max_diff:
            max_diff = np.abs(conv_diff).max()

        if len(conv_diff) <= 9: 
            conv_diff = diff
        
        bins = bins - 0.5
        bins = bins[1:]
        scatter_data = pd.DataFrame({
            'bin_midpoints': bins,
            'diff': conv_diff
        })
        scatter_data = scatter_data.loc[scatter_data['diff'] != 0]

        data_dict[loc] = (df_melt, scatter_data)

    return data_dict, max_diff

def plot_data(data_dict, max_diff, num_cols, cancer, motifs, sample_case, sample_control):
    num_loci = len(data_dict.keys())
    num_rows = np.ceil(num_loci / num_cols).astype(int)
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(24, 3+num_rows*5))
    axs = axs.ravel()

    plt.subplots_adjust(wspace=0.3, hspace=0.5)


    for i, (loc, (df_melt, scatter_data)) in enumerate(data_dict.items()):
        sns.histplot(data=df_melt, x='Repeat Length', hue='group', element='step', common_norm=False, binwidth=1, ax=axs[i])
    
        ax2 = axs[i].twinx()  # Create a second axes that shares the same x-axis
        ax2.set_ylim(-max_diff * 2, max_diff * 2)

        sns.lineplot(data=scatter_data, x='bin_midpoints', y='diff', ax=ax2, color='purple')
        ax2.axhline(0, color='gray', linestyle='--')
        ax2.set_ylabel('Smoothed Case - Control Difference')
        ax2.set_xlim(0, scatter_data['bin_midpoints'].max() + 2)

        legend_elements = [Patch(facecolor='b', edgecolor='b', label='Control', alpha=0.3),
                            Patch(facecolor='orange', edgecolor='orange', label='Case', alpha=0.3),
                            Line2D([0], [0], color='purple', label='Smoothed Difference'),
                            Line2D([0], [0], color='gray', linestyle='--', label='Zero Difference')]

        if sample_case is not None:
            vals = sample_case.loc[:, loc].dropna().values
            for val in vals:
                axs[i].axvline(val, color='red', linestyle='--')
            legend_elements.append(Line2D([0], [0], color='red', linestyle='--', label='Sample Case Length'))
            

        if sample_control is not None:
            vals = sample_control.loc[:, loc].dropna().values
            for val in vals:
                axs[i].axvline(val, color='blue', linestyle='--')
            legend_elements.append(Line2D([0], [0], color='blue', linestyle='--', label='Sample Control Length'))

        axs[i].legend(handles=legend_elements, loc='upper right')
        axs[i].set_title(f'{loc}, {motifs["chr" + loc]}')
        # add vertical line at sample repeat length
        

        axs[i].set_xlabel('Repeat Length')

    plt.suptitle(f'Genotypes of {cancer} Across Loci')
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.show()

def graphLociGenotypes(loci, case_df, control_df, cancer, sample=None, num_cols=4):
    motifs = get_motifs(loci)

    sample_cases = case_df[case_df['sample_id'] == sample] if sample else None
    sample_controls = control_df[control_df['sample_id'] == sample] if sample else None

    data_dict, max_diff = calculate_data(loci, case_df, control_df)
    plot_data(data_dict, max_diff, num_cols, cancer, motifs, sample_cases, sample_controls)

def get_motifs(loci):
    df = pd.read_csv('data/other/locus_structures.csv')

    # add 'chr' to start of loci
    loci = ['chr' + loc for loc in loci]

    # return dict of loci in loci and their motifs
    df = df.loc[df['ReferenceRegion'].isin(loci)]
    return dict(zip(df['ReferenceRegion'], df['LocusStructure']))



def graphCancersGenotypes(locus, cancers):
    """
    Plots the genotypes of a single locus for each cancer type.
    args:
        locus: the locus to plot
        cancers: a list of cancer types to plot
    """
    num_cancers = len(cancers)
    num_cols = 4
    num_rows = np.ceil(num_cancers / num_cols).astype(int)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, num_rows*5))  # Adjust the figure size as per your requirement
    axs = axs.ravel()  # flatten the array of axes

    # Adjust padding between subplots
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.subplots_adjust(top=0.95)  # adjust top margin

    for i, cancer in enumerate(cancers):
        # Read in the case and control data for each cancer
        case_df = pd.read_csv(f'data/genotypes/{cancer}_case_tidy.csv', usecols=[locus], index_col=False)
        control_df = pd.read_csv(f'data/genotypes/{cancer}_control_tidy.csv', usecols=[locus], index_col=False)

        # Combine the case and control dataframes along columns
        df = pd.concat([case_df, control_df], axis=1)
        df.columns = ['case', 'control']
        df.dropna(inplace=True)

        # Melt the dataframe to long format for plotting
        df_melt = df.melt(value_name='Repeat Length', var_name='group')
        bins = np.arange(df_melt['Repeat Length'].min() - 0.5, df_melt['Repeat Length'].max() + 0.5)
        # Create the overlaid histograms
        sns.histplot(data=df_melt, x='Repeat Length', hue='group', element='step', common_norm=False, bins=bins, ax=axs[i])

        
        # Get histogram data for case and control groups
        hist_case, bins_case = np.histogram(df['case'], bins=bins)
        hist_control, bins_control = np.histogram(df['control'], bins=bins)

        # Calculate the difference in bin counts
        diff = hist_case - hist_control

        max_hist = np.max([hist_case.max(), hist_control.max()])
        ax2 = axs[i].twinx()  # Create a second axes that shares the same x-axis
        ax2.set_ylim(-max_hist/2, max_hist/2)  # Set the limits of the y-axis to the maximum difference

        # Prepare data for the smoothed line plot
        scatter_data = pd.DataFrame({
            'bin_midpoints': (bins_case[:-1] + bins_case[1:])/2,
            'diff': diff
        })
        scatter_data = scatter_data.loc[scatter_data['diff'] != 0]

        sns.regplot(data=scatter_data, x='bin_midpoints', y='diff', order=4, ax=ax2, color='purple', scatter=False, ci = None, truncate=True)
        ax2.set_ylabel('Difference in Bin Counts')
        ax2.set_xlim(0, scatter_data['bin_midpoints'].max() + 2)

        axs[i].set_title(locus)
        axs[i].set_xlabel('Repeat Length')

    plt.suptitle(f'Genotype Distributions of {locus} across Cancers')
    plt.tight_layout()
    plt.show()





### DIFFERENCE GRAPHING FUNCTIONS ###

def graphDiff(loc, cancer):
    dat = pd.read_csv(f'data/diffs/{cancer}_diff.csv', usecols=[loc])
    dat[loc].hist(bins=100)
    plt.title(f'{cancer} ; {loc} Distribution')
    plt.xlabel('Tumor - Normal')
    plt.ylabel('Frequency')
    # Comment out not to save the figure
    # plt.savefig(f'figs/{cancer}_{loc}.png')
    plt.show()

def graphDiffLociHelper(loc, cancer, ax):
    dat = pd.read_csv(f'data/diffs/{cancer}_diff.csv', usecols=[loc])
    ax.hist(dat[loc], bins=100)
    ax.set_title(f'{loc} Distribution')
    ax.set_xlabel('Tumor - Normal')
    ax.set_ylabel('Frequency')
    

def graphLociDiffs(loci, cancer):
    num_loci = len(loci)
    num_cols = 4
    num_rows = math.ceil(num_loci / num_cols)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(30, num_rows*5)) # Adjust the figure size as per your requirement
    axs = axs.ravel() # flatten the array of axes

    for idx, loc in enumerate(loci):
        ax = axs[idx]
        graphDiffLociHelper(loc, cancer, ax)

    plt.title(f'{cancer} Differences (Tumor-Normal) Distributions')
    plt.tight_layout()  # To prevent overlap of subplots
    plt.show()






def graphLoci(loc, cancers):
    # Assuming cancers is a list of your cancers
    num_cancers = len(cancers)
    num_cols = 5
    num_rows = math.ceil(num_cancers / num_cols)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(30, num_rows*5)) # Adjust the figure size as per your requirement
    axs = axs.ravel() # flatten the array of axes

    for idx, cancer in enumerate(cancers):
        ax = axs[idx]
        try:
            graphLocus(loc, cancer, ax)
        except:
            print(f'Could not graph {cancer}')
            continue

    plt.tight_layout()  # To prevent overlap of subplots
    plt.show()

def graphLocus(loc, cancer, ax):
    dat = pd.read_csv(f'data/diffs/{cancer}_diff.csv', usecols=[loc])
    ax.hist(dat[loc], bins=100)
    ax.set_title(f'{cancer}Distribution')
    ax.set_xlabel('Tumor - Normal')
    ax.set_ylabel('Frequency')

