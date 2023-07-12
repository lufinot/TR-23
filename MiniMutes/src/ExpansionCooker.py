import argparse
import pandas as pd
import orjson
import os
import numpy as np
import logging

MAX_WIDTH = 2

def is_wide(ci):
    ci = list(map(int, ci.split('-')))
    return ci[1] - ci[0] > MAX_WIDTH

def decide_genotype_order(case, control):
    if case[1] == control[0]:
        temp = case[1]
        case[1] = case[0]
        case[0] = temp


    return case, control

def extract_genotypes_diffs(manifest_path, disease_name, raw_eh_dir, output_dir):

    # Load the manifest file
    manifest = pd.read_csv(manifest_path)

    # Ensure output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Create a log file
    logging.basicConfig(filename=os.path.join(output_dir, f'{disease_name}_genotype_getter.log'), level=logging.INFO)

    # initalize empty dataframes for case and control
    genotype_df = []
    df_tracking = []
    for i, row in manifest.iterrows():
        icgc_donor_id = row['icgc_donor_id']
        file_path_case = os.path.join(raw_eh_dir, f"{row['case_object_id']}.json")
        file_path_control = os.path.join(raw_eh_dir, f"{row['control_object_id']}.json")  
        with open(file_path_case, 'r') as file_case, open(file_path_control, 'r') as file_control:
            try:
                data_case = orjson.loads(file_case.read())
                data_control = orjson.loads(file_control.read())
            except:
                logging.error(f"Error loading {icgc_donor_id}")
                continue
            for locus in set(data_case['LocusResults']):
                allele_count = data_case['LocusResults'][locus]['AlleleCount']
                for variant in set(data_case['LocusResults'][locus]['Variants']):
                    case = data_case['LocusResults'][locus]['Variants'][variant]
                    control = data_control['LocusResults'][locus]['Variants'][variant]
                    refRegion = data_case['LocusResults'][locus]['Variants'][variant]['ReferenceRegion']

                    # get genotypes
                    case_ci = case.get('GenotypeConfidenceInterval')
                    control_ci = control.get('GenotypeConfidenceInterval')
                    if case_ci is None or control_ci is None:
                        continue

                    if allele_count == 1:
                        if is_wide(case_ci) or is_wide(control_ci):
                            df_tracking.append({'icgc_donor_id': icgc_donor_id, 
                                                'locus': locus, 
                                                'motif': case.get('RepeatUnit'),
                                                'control_ci': control_ci, 
                                                'case_ci': case_ci})
                            continue
                        control_genotype = int(control.get('Genotype'))
                        case_genotype = int(case.get('Genotype'))
                        genotype_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'motif': case.get('RepeatUnit'), 'control': control_genotype , 'case': case_genotype, 'diff': case_genotype - control_genotype})
                        
                        continue
                    if allele_count == 2:
                        case_ci = case_ci.split('/')
                        control_ci = control_ci.split('/')
                        case_genotypes = case.get('Genotype')
                        control_genotypes = control.get('Genotype')
                        case_genotypes = list(map(int, case_genotypes.split('/')))
                        control_genotypes = list(map(int, control_genotypes.split('/')))
                        
                        if is_wide(case_ci[0]) or is_wide(control_ci[0]):
                            df_tracking.append({'icgc_donor_id': icgc_donor_id, 
                                                'locus': locus, 
                                                'motif': case.get('RepeatUnit'),
                                                'control_ci': control_ci[0], 
                                                'case_ci': case_ci[0]})
                            case_genotypes[0] = np.nan
                            control_genotypes[0] = np.nan
                        
                        if is_wide(case_ci[1]) or is_wide(control_ci[1]):
                            df_tracking.append({'icgc_donor_id': icgc_donor_id, 
                                                'locus': locus, 
                                                'motif': case.get('RepeatUnit'),
                                                'control_ci': control_ci[1], 
                                                'case_ci': case_ci[1]})
                            case_genotypes[1] = np.nan
                            control_genotypes[1] = np.nan

                        if (case_genotypes[0] is None and case_genotypes[1] is None) or (control_genotypes[0] is None and control_genotypes[1] is None):
                            continue

                        control_genotypes, case_genotypes = decide_genotype_order(control_genotypes, case_genotypes)
                        genotype_df.append({'donor_id': icgc_donor_id + '_0', 'refRegion': refRegion, 'motif': case.get('RepeatUnit'),'control': control_genotypes[0], 'case': case_genotypes[0], 'diff': case_genotypes[0] - control_genotypes[0]})
                        genotype_df.append({'donor_id': icgc_donor_id + '_1', 'refRegion': refRegion, 'motif': case.get('RepeatUnit'),'control': control_genotypes[1], 'case': case_genotypes[1], 'diff': case_genotypes[1] - control_genotypes[1]})

    genotype_df = pd.DataFrame(genotype_df)
    df_tracking = pd.DataFrame(df_tracking)


    # Save the output dataframes
    genotype_df.to_csv(os.path.join(output_dir, f'{disease_name}_genotypes_diffs.csv'))
    df_tracking.to_csv(os.path.join(output_dir, f'{disease_name}_tracking.csv'))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process BAMlet files.')
    parser.add_argument('--manifest', required=True, help='Path to the manifest file.')
    parser.add_argument('--disease', required=True, help='Disease name for the run of the script.')
    parser.add_argument('--raw_eh', required=True, help='Directory with raw EH.')
    parser.add_argument('--output', required=True, help='Directory name for the output (created if not existing).')

    args = parser.parse_args()

    extract_genotypes_diffs(args.manifest, args.disease, args.raw_eh, args.output)
