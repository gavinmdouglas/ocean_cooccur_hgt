#!/usr/bin/python3

import os
import sys
from collections import defaultdict

# Parse folders with FASTP-filtered FASTQs and:
# (1) For run IDs that are the sole runs for a sample, produce symolic links to renamed FASTQs that include the sample ID instead.
# (2) For run IDs that are one of multiple of the same sample, produce cat commands to concatenate the FASTQ(s) for that sample.

nico_dir = dict()
nico_dir['PRJEB1787'] = '/mfs/nicot/PRJEB1787-fastq_ftp-20231101-1507/cleaned/'
nico_dir['PRJEB9740'] = '/mfs/nicot/PRJEB9740-fastq_ftp-20231101-1513/cleaned/'

# Read through Tara ocean metadata tables first.
# This RNA-seq table needs to be run separately: 'PRJEB6608'
TARA_bioprojects = ['PRJEB1787', 'PRJEB9740']
TARA_sample_to_run = defaultdict(list)
sample_to_bioproject = dict()
for TARA_bioproject in TARA_bioprojects:
    metadata_file = '/mfs/gdouglas/projects/ocean_mags/metadata/' + TARA_bioproject + '_metadata.csv'
    with open(metadata_file, 'r') as metadata_fh:
        metadata_header = metadata_fh.readline().strip().split(',')
        metadata_col_to_i = {col: i for i, col in enumerate(metadata_header)}
        run_col_i = metadata_col_to_i['Run']
        sample_col_i = metadata_col_to_i['Sample']
        for line in metadata_fh:
            line = line.strip().split(',')
            run_ID = line[run_col_i]

            # Skip run IDs that are from same sample for which there are already technical replicates that are PE
            # (easier just to stick with the PE samples).
            if run_ID in set(['ERR164409', 'ERR164408', 'ERR164407']):
                continue

            sample_ID = line[sample_col_i]
            TARA_sample_to_run[sample_ID].append(run_ID)

            if sample_ID in sample_to_bioproject.keys():
                if sample_to_bioproject[sample_ID] != TARA_bioproject:
                    sys.exit("Error: sample ID in multiple bioprojects: " + sample_ID)
            else:
                sample_to_bioproject[sample_ID] = TARA_bioproject

num_unique_reps = 0
num_multi_reps = 0
ln_cmds = []
cat_cmds = []

for sample_id in sample_to_bioproject.keys():
    runs = TARA_sample_to_run[sample_id]
    if len(runs) == 0:
        sys.exit("Error: no runs for " + sample_id + " ?")
    elif len(runs) == 1:
        num_unique_reps += 1
        RUNID = runs[0]
        exp_R1_try1 = nico_dir[sample_to_bioproject[sample_id]] + RUNID + '_1_clean.fastq.gz'
        exp_R2_try1 = nico_dir[sample_to_bioproject[sample_id]] + RUNID + '_2_clean.fastq.gz'

        exp_R1_try2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE/' + RUNID + '_1.fastq.gz'
        exp_R2_try2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE/' + RUNID + '_2.fastq.gz'

        out_R1 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_Tara/' + sample_id + '_1.fastq.gz'
        out_R2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_Tara/' + sample_id + '_2.fastq.gz'

        if os.path.exists(exp_R1_try1) and os.path.exists(exp_R2_try1) and os.path.exists(exp_R1_try2) and os.path.exists(exp_R2_try2):
            sys.exit("Error: FASTQ file in both places " + exp_R1_try1)
        elif os.path.exists(exp_R1_try1) and os.path.exists(exp_R2_try1):
            ln_cmds.append("ln -s " + exp_R1_try1 + " " + out_R1)
            ln_cmds.append("ln -s " + exp_R2_try1 + " " + out_R2)
        elif os.path.exists(exp_R1_try2) and os.path.exists(exp_R2_try2):
            ln_cmds.append("ln -s " + exp_R1_try2 + " " + out_R1)
            ln_cmds.append("ln -s " + exp_R2_try2 + " " + out_R2)
        else:
            sys.exit("Error: FASTQ file not found " + exp_R1_try1)
    elif len(runs) > 1:
        num_multi_reps += 1

        cat_cmd_R1 = "cat "
        cat_cmd_R2 = "cat "
        cat_cmd_SE = "cat "
        pe_counts = 0
        se_counts = 0
        for RUNID in runs:
            exp_R1_try1 = nico_dir[sample_to_bioproject[sample_id]] + RUNID + '_1_clean.fastq.gz'
            exp_R2_try1 = nico_dir[sample_to_bioproject[sample_id]] + RUNID + '_2_clean.fastq.gz'

            exp_R1_try2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE/' + RUNID + '_1.fastq.gz'
            exp_R2_try2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE/' + RUNID + '_2.fastq.gz'

            if os.path.exists(exp_R1_try1) and os.path.exists(exp_R2_try1) and os.path.exists(exp_R1_try2) and os.path.exists(exp_R2_try2):
                sys.exit("Error: FASTQ file in both places " + exp_R1_try1)
            elif os.path.exists(exp_R1_try1) and os.path.exists(exp_R2_try1):
                cat_cmd_R1 += " " + exp_R1_try1
                cat_cmd_R2 += " " + exp_R2_try1
                pe_counts += 1
            elif os.path.exists(exp_R1_try2) and os.path.exists(exp_R2_try2):
                cat_cmd_R1 += " " + exp_R1_try2
                cat_cmd_R2 += " " + exp_R2_try2
                pe_counts += 1
            else:
                # Try SE.
                exp_SE_try1 = nico_dir[sample_to_bioproject[sample_id]] + RUNID + '_clean.fastq.gz'
                exp_SE_try2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_SE/' + RUNID + '.fastq.gz'
                if os.path.exists(exp_SE_try1) and os.path.exists(exp_SE_try2):
                    sys.exit("Error: FASTQ file in both places " + exp_SE_try1)
                elif os.path.exists(exp_SE_try1) :
                    cat_cmd_SE += " " + exp_SE_try1
                    se_counts += 1
                elif os.path.exists(exp_SE_try2):
                    cat_cmd_SE += " " + exp_SE_try2
                    se_counts += 1
        
        if pe_counts > 0 and se_counts > 0:
            print(cat_cmd_SE, file=sys.stderr)
            print(cat_cmd_R1, file=sys.stderr)
            print(cat_cmd_R2, file=sys.stderr)
            print(runs, file=sys.stderr)
            sys.exit("Error: both PE and SE files found for " + sample_id)
        elif pe_counts > 0 and se_counts == 0:
            out_R1 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_Tara/' + sample_id + '_1.fastq.gz'
            out_R2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_Tara/' + sample_id + '_2.fastq.gz'
            cat_cmd_R1 += " > " + out_R1
            cat_cmd_R2 += " > " + out_R2
            cat_cmds.append(cat_cmd_R1)
            cat_cmds.append(cat_cmd_R2)
        elif pe_counts == 0 and se_counts > 0:
            out_SE = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_SE_final_Tara/' + sample_id + '.fastq.gz'
            cat_cmd_SE += " > " + out_SE
            cat_cmds.append(cat_cmd_SE)
        else:
            sys.exit("Error: no PE or SE files found for " + sample_id)

print("Unique reps:", file=sys.stderr)
print(num_unique_reps, file=sys.stderr)

print("Multi reps:", file=sys.stderr)
print(num_multi_reps, file=sys.stderr)

# Write out commands.
with open("/mfs/gdouglas/projects/ocean_mags/additional_round2/Tara_ln_cmds.sh", "w") as ln_fh:
    for ln_cmd in ln_cmds:
        print(ln_cmd, file=ln_fh)

with open("/mfs/gdouglas/projects/ocean_mags/additional_round2/Tara_cat_cmds.sh", "w") as cat_fh:
    for cat_cmd in cat_cmds:
        print(cat_cmd, file=cat_fh)
