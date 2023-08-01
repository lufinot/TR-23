import argparse
import pandas as pd
import orjson
import os
import numpy as np
import logging

from multiprocessing import Pool, cpu_count
from functools import partial

MAX_WIDTH = 5

LOG_LEVEL = os.getenv('LOG_LEVEL') or 'info'

log_dict = {'debug':logging.DEBUG, 'info':logging.INFO , 'warning':logging.WARNING, 
            'error':logging.ERROR, 'critical':logging.CRITICAL}

log_level = log_dict.get(LOG_LEVEL.lower(),logging.INFO)

logging.basicConfig(filename='filename.log', level=log_level,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

def is_wide(ci):
    ci = list(map(int, ci.split('-')))
    return ci[1] - ci[0] > MAX_WIDTH

def decide_genotype_order(case, control):

    if case[1] == control[0]:
        temp = case[1]
        case[1] = case[0]
        case[0] = temp


    return case, control

def getPairs(case_ci, control_ci, case_genotypes, control_genotypes):
    if is_wide(case_ci[0]):
        if is_wide(control_ci[1]):
            return [case_genotypes[1], control_genotypes[0]], [case_ci[0], control_ci[1]]
        else:
            return [case_genotypes[1], control_genotypes[1]], [case_ci[0], control_ci[0]]
    else:
        if is_wide(control_ci[0]):
            return [case_genotypes[0], control_genotypes[1]], [case_ci[1], control_ci[0]]
        else:
            return [case_genotypes[0], control_genotypes[0]], [case_ci[1], control_ci[1]]
        
def process_locus(icgc_donor_id, data_case, data_control, local_case_df, local_control_df, local_diff_df, local_df_tracking):
    allele_count = data_case['AlleleCount']
    for variant in set(data_case['Variants']):
        case = data_case['Variants'][variant]
        control = data_control['Variants'][variant]
        refRegion = data_case['Variants'][variant]['ReferenceRegion']

        # get values
        case_ci = case.get('GenotypeConfidenceInterval')
        control_ci = control.get('GenotypeConfidenceInterval')
        if case_ci is None or control_ci is None:
            continue

        if allele_count == 1:
            if is_wide(case_ci) or is_wide(control_ci):
                local_df_tracking.append({'icgc_donor_id': icgc_donor_id, 
                                    'refRegion': refRegion,
                                    'motif': case.get('RepeatUnit'),
                                    'control_ci': control_ci, 
                                    'case_ci': case_ci})
                continue

            local_case_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': case.get('Genotype')})
            local_control_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': control.get('Genotype')})
            continue
        if allele_count == 2:
            case_ci = case_ci.split('/')
            control_ci = control_ci.split('/')
            case_genotypes = case.get('Genotype')
            control_genotypes = control.get('Genotype')
            case_genotypes = list(map(int, case_genotypes.split('/')))
            control_genotypes = list(map(int, control_genotypes.split('/')))
            # make any genotypes witha  wide confidence interval nan
            if is_wide(case_ci[0]):
                case_genotypes[0] = np.nan
            if is_wide(case_ci[1]):
                case_genotypes[1] = np.nan
            if is_wide(control_ci[0]):
                control_genotypes[0] = np.nan
            if is_wide(control_ci[1]):
                control_genotypes[1] = np.nan

            
            tot_wide = np.isnan(case_genotypes).sum() + np.isnan(control_genotypes).sum()

            if tot_wide == 0:
                case_genotypes, control_genotypes = decide_genotype_order(case_genotypes, control_genotypes)
                local_case_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': case_genotypes[0]})
                local_case_df.append({'donor_id': icgc_donor_id + '_1', 'refRegion': refRegion, 'value': case_genotypes[1]})
                local_control_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': control_genotypes[0]})
                local_control_df.append({'donor_id': icgc_donor_id + '_1', 'refRegion': refRegion, 'value': control_genotypes[1]})
                local_diff_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': case_genotypes[0] - control_genotypes[0]})
                local_diff_df.append({'donor_id': icgc_donor_id + '_1', 'refRegion': refRegion, 'value': case_genotypes[1] - control_genotypes[1]})
                continue

            # if 3 out of 4 values are nan, skip
            if tot_wide >= 3 or (is_wide(case_ci[0]) and is_wide(case_ci[1])) or (is_wide(control_ci[0]) and is_wide(control_ci[1])):
                local_df_tracking.append({'icgc_donor_id': icgc_donor_id, 
                                    'refRegion': refRegion, 
                                    'motif': case.get('RepeatUnit'),
                                    'control_ci': control_ci[0], 
                                    'case_ci': case_ci[0]})
                local_df_tracking.append({'icgc_donor_id': icgc_donor_id, 
                                    'refRegion': refRegion, 
                                    'motif': case.get('RepeatUnit'),
                                    'control_ci': control_ci[1], 
                                    'case_ci': case_ci[1]})
                continue
            
            goodPair, badCi = getPairs(case_ci, control_ci, case_genotypes, control_genotypes)
            local_case_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': goodPair[0]})
            local_control_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': goodPair[1]})
            local_diff_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'value': goodPair[0] - goodPair[1]})
            local_df_tracking.append({'icgc_donor_id': icgc_donor_id, 
                                'refRegion': refRegion,
                                'motif': case.get('RepeatUnit'),
                                'control_ci': badCi[1], 
                                'case_ci': badCi[0]})


def process_donor(donor, raw_eh_dir):
    icgc_donor_id = donor['icgc_donor_id']
    file_path_case = os.path.join(raw_eh_dir, f"{donor['case_object_id']}.json")
    file_path_control = os.path.join(raw_eh_dir, f"{donor['control_object_id']}.json")

    local_case_df = []
    local_control_df = []
    local_diff_df = []
    local_df_tracking = []

    # test if the files exist
    if not os.path.isfile(file_path_case) or not os.path.isfile(file_path_control):
        logging.error(f'One or more files do not exist: {file_path_case}, {file_path_control}')
        return local_case_df, local_control_df, local_diff_df, local_df_tracking, icgc_donor_id

    with open(file_path_case, 'r') as file_case, open(file_path_control, 'r') as file_control:
        try:
            data_case = orjson.loads(file_case.read())
            data_control = orjson.loads(file_control.read())
        except Exception as e:
            logging.error(f'Error decoding JSON for files {file_path_case}, {file_path_control}. Error: {str(e)}')
            return local_case_df, local_control_df, local_diff_df, local_df_tracking

        for locus in set(data_case['LocusResults']):
            process_locus(icgc_donor_id, data_case['LocusResults'][locus], data_control['LocusResults'][locus], local_case_df, local_control_df, local_diff_df, local_df_tracking)


    return local_case_df, local_control_df, local_diff_df, local_df_tracking, icgc_donor_id


def extract_genotypes_diffs(manifest_path, disease_name, raw_eh_dir, output_dir):
    # Load the manifest
    manifest = pd.read_csv(manifest_path)

    # Ensure output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)


    # Use multiprocessing to process each donor's files
    with Pool(processes=cpu_count()) as pool:
        func = partial(process_donor, raw_eh_dir=raw_eh_dir)
        results = pool.map(func, manifest.to_dict('records'))

    # Aggregate results
    case_df_list, control_df_list, diff_df_list, df_tracking_list, donor_ids = zip(*results)
    case_df = [item for sublist in case_df_list for item in sublist]
    control_df = [item for sublist in control_df_list for item in sublist]
    diff_df = [item for sublist in diff_df_list for item in sublist]
    df_tracking = [item for sublist in df_tracking_list for item in sublist]

    # Convert lists to DataFrames and pivot
    case_df = pd.DataFrame(case_df).pivot(index='donor_id', columns='refRegion', values='value')
    control_df = pd.DataFrame(control_df).pivot(index='donor_id', columns='refRegion', values='value')
    diff_df = pd.DataFrame(diff_df).pivot(index='donor_id', columns='refRegion', values='value')
    df_tracking = pd.DataFrame(df_tracking)

    # Save the DataFrames
    case_df.to_csv(os.path.join(output_dir, f'{disease_name}_case.csv'))
    control_df.to_csv(os.path.join(output_dir, f'{disease_name}_control.csv'))
    diff_df.to_csv(os.path.join(output_dir, f'{disease_name}_diff.csv'))
    df_tracking.to_csv(os.path.join(output_dir, f'{disease_name}_tracking.csv'))

    return case_df, control_df, diff_df, df_tracking


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process BAMlet files.')
    parser.add_argument('--manifest', required=True, help='Path to the manifest file.')
    parser.add_argument('--disease', required=True, help='Disease name for the run of the script.')
    parser.add_argument('--raw_eh', required=True, help='Directory with raw EH.')
    parser.add_argument('--output', required=True, help='Directory name for the output (created if not existing).')

    args = parser.parse_args()

    _, _, _, _ = extract_genotypes_diffs(args.manifest, args.disease, args.raw_eh, args.output)
