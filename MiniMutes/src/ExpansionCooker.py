import argparse
import pandas as pd
import orjson
import os
import numpy as np
import logging
from datetime import datetime
import logging.handlers
import queue

from multiprocessing import Pool, cpu_count
from functools import partial

MAX_WIDTH = 5




LOG_LEVEL = os.getenv('LOG_LEVEL') or 'info'
log_dict = {'debug': logging.DEBUG, 'info': logging.INFO, 'warning': logging.WARNING, 
            'error': logging.ERROR, 'critical': logging.CRITICAL}
log_level = log_dict.get(LOG_LEVEL.lower(), logging.INFO)

slurm_job_id = os.getenv('SLURM_JOBID')
if not slurm_job_id:
    slurm_job_id = datetime.now().strftime('%Y%m%d_%H%M')


log_dir = 'ExpansionCookerLogs'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

logging.basicConfig(filename=os.path.join(log_dir, f'{slurm_job_id}_ExpansionCooker.log'), 
                    level=log_level,
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

class GenotypeChecker:
    def __init__(self, genotypes, spanning_reads, flanking_reads):
        self.genotypes = list(map(int, genotypes.split('/')))
        self.spanning_reads_dict = self._parse_reads_to_dict(spanning_reads)
        self.flanking_reads_list = sorted(self.parse_reads_to_list(flanking_reads))
        self.supported_genotypes = []


def process_locus(donor_id, data_case, data_control, local_case_df, local_control_df, local_diff_df, local_df_tracking):
    for variant in set(data_case['Variants']):
        case = data_case['Variants'][variant]
        control = data_control['Variants'][variant]
        refRegion = data_case['Variants'][variant]['ReferenceRegion']

        control_genotypes = list(map(int, control.get('Genotype').split('/')))
        case_genotypes = list(map(int, case.get('Genotype').split('/')))
        
        is_control_homozygoous = len(control_genotypes) == 2 and control_genotypes[0] == control_genotypes[1]

        if is_control_homozygoous:
            # treat it as a single allele
            control_genotypes = [control_genotypes[0]]
      
        # make genotype checker objects for case and cotrol
        case_genotype_checker = GenotypeChecker(case_genotypes, case.get('SpanningReads'), case.get('FlankingReads'))
        control_genotype_checker = GenotypeChecker(control_genotypes, control.get('SpanningReads'), control.get('FlankingReads'))

        # get supported genotypes for control

        # if control is homozygous, add the same genotype twice to the supported genotypes

        # add control's supported genotypes to beggining of case genotypes

        # get supported genotypes for case, with the list that now has control's supported genotypes

        
        # match the genotypes 

        # add to lists

def process_locus(donor_id, data_case, data_control, local_case_df, local_control_df, local_diff_df, local_df_tracking):
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
                local_df_tracking.append({'donor_id': donor_id, 
                                    'refRegion': refRegion,
                                    'motif': case.get('RepeatUnit'),
                                    'control_ci': control_ci, 
                                    'case_ci': case_ci})
                continue

            local_case_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': case.get('Genotype')})
            local_control_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': control.get('Genotype')})
            continue
        if allele_count == 2:

            # check number of supporting reads for case and control
            # if either has a low number (gotta choose threshold) 
            #   add to tracking or  
            # If theres a big difference in amounts
            #   Figure out rubost way to deal with this
            # Else
            #   Trust the values

            case_ci = case_ci.split('/')
            control_ci = control_ci.split('/')
            case_genotypes = case.get('Genotype')
            control_genotypes = control.get('Genotype')
            case_genotypes = list(map(int, case_genotypes.split('/')))
            control_genotypes = list(map(int, control_genotypes.split('/')))
            # make any genotypes with a wide confidence interval nan
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
                local_case_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': case_genotypes[0]})
                local_case_df.append({'donor_id': donor_id + '_1', 'refRegion': refRegion, 'value': case_genotypes[1]})
                local_control_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': control_genotypes[0]})
                local_control_df.append({'donor_id': donor_id + '_1', 'refRegion': refRegion, 'value': control_genotypes[1]})
                local_diff_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': case_genotypes[0] - control_genotypes[0]})
                local_diff_df.append({'donor_id': donor_id + '_1', 'refRegion': refRegion, 'value': case_genotypes[1] - control_genotypes[1]})
                continue

            # if 3 out of 4 values are nan, skip
            if tot_wide >= 3 or (is_wide(case_ci[0]) and is_wide(case_ci[1])) or (is_wide(control_ci[0]) and is_wide(control_ci[1])):
                local_df_tracking.append({'donor_id': donor_id, 
                                    'refRegion': refRegion, 
                                    'motif': case.get('RepeatUnit'),
                                    'control_ci': control_ci[0], 
                                    'case_ci': case_ci[0]})
                local_df_tracking.append({'donor_id': donor_id, 
                                    'refRegion': refRegion, 
                                    'motif': case.get('RepeatUnit'),
                                    'control_ci': control_ci[1], 
                                    'case_ci': case_ci[1]})
                continue
            
            goodPair, badCi = getPairs(case_ci, control_ci, case_genotypes, control_genotypes)
            local_case_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': goodPair[0]})
            local_control_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': goodPair[1]})
            local_diff_df.append({'donor_id': donor_id + '_0', 'refRegion': refRegion, 'value': goodPair[0] - goodPair[1]})
            local_df_tracking.append({'donor_id': donor_id, 
                                'refRegion': refRegion,
                                'motif': case.get('RepeatUnit'),
                                'control_ci': badCi[1], 
                                'case_ci': badCi[0]})


def process_donor(donor, raw_eh_dir):
    donor_id = donor['donor_id']
    logging.info(f'Processing {donor_id}.')
    file_path_case = os.path.join(raw_eh_dir, f"{donor['case_object_id']}.json")
    file_path_control = os.path.join(raw_eh_dir, f"{donor['control_object_id']}.json")

    local_case_df = []
    local_control_df = []
    local_diff_df = []
    local_df_tracking = []

       # Test if the files exist
    case_exists = os.path.isfile(file_path_case)
    control_exists = os.path.isfile(file_path_control)

    if not case_exists and not control_exists:
        logging.error(f'Missing both files for{donor_id}, both files do not exist: {file_path_case}, {file_path_control}')
    elif not case_exists:
        logging.error(f'Missing case for {donor_id}: {file_path_case}')
    elif not control_exists:
        logging.error(f'Missing control for {donor_id}: {file_path_control}')
        
    if not case_exists or not control_exists:
        return local_case_df, local_control_df, local_diff_df, local_df_tracking, donor_id

    with open(file_path_case, 'r') as file_case, open(file_path_control, 'r') as file_control:
        try:
            data_case = orjson.loads(file_case.read())
            data_control = orjson.loads(file_control.read())
        except Exception as e:
            logging.error(f'Could not decode JSON for {donor_id} Error: {str(e)}')
            return local_case_df, local_control_df, local_diff_df, local_df_tracking, donor_id

        for locus in set(data_case['LocusResults']):
            process_locus(donor_id, data_case['LocusResults'][locus], data_control['LocusResults'][locus], local_case_df, local_control_df, local_diff_df, local_df_tracking)

    logging.info(f'Finished {donor_id} succesfully.')
    return local_case_df, local_control_df, local_diff_df, local_df_tracking, donor_id


def extract_genotypes_diffs(manifest_path, disease_name, raw_eh_dir, output_dir):
    
    # Load the manifest
    manifest = pd.read_csv(manifest_path)


    # Ensure output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    with Pool(processes=cpu_count()) as pool:
        func = partial(process_donor, raw_eh_dir=raw_eh_dir)
        results = pool.map(func, manifest.to_dict('records'))


    logging.info('Finished processing files, combining results.')
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

    logging.info(f'Saving DataFrames.')

    # Save the DataFrames
    case_df.to_csv(os.path.join(output_dir, f'{disease_name}_case.csv'))
    control_df.to_csv(os.path.join(output_dir, f'{disease_name}_control.csv'))
    diff_df.to_csv(os.path.join(output_dir, f'{disease_name}_diff.csv'))
    df_tracking.to_csv(os.path.join(output_dir, f'{disease_name}_tracking.csv'))

    logging.info('Finished Saving DataFrames.')

    return diff_df

def init_argparse():
    parser = argparse.ArgumentParser(description='Process Expansion Hunter output for analysis of paired genotype differences.')
    parser.add_argument('raw_eh', metavar='RawDir', type=str, help='Directory with Expansion Hunter output JSONs.')
    parser.add_argument('manifest', metavar='Manifest', type=str, help='Manifest file with case and control object ids.')
    parser.add_argument('--name', '-n', required=True, help='Disease name for output files.')
    parser.add_argument('--outdir', '-o', required=True, help='Output directory (default .).')
    parser.add_argument('--feats', '-f', default=False, action='store_true', help='Create features from the output? (Default: False)')
    return parser


def main(args=None):
    parser = init_argparse()
    args = parser.parse_args()

    diffs = extract_genotypes_diffs(args.manifest, args.name, args.raw_eh, args.outdir)

    if args.feats:
        logging.info('Creating features from the output.')
        # import EH_Feature_Extractor as efe
        # feats = process_and_extract_features(diffs)
        # feats.to_csv(os.path.join(args.out, f'{args.name}_feats.csv'), index=False)
        logging.info('Finished creating features.')


if __name__ == "__main__":
    main()
    
  