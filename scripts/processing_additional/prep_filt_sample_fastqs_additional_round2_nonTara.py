#!/usr/bin/python3

import os
import sys
from collections import defaultdict

# Read through OceanDNA samples excluded last time due to being from same biological sample, but different filter cut-off
# (I want to process everything this time).

sample_to_runs = dict()
tmp_runs = set()
sample_to_folder = dict()
sample_to_seqtype = dict()

pe_runs = set()
with open('/mfs/gdouglas/projects/ocean_mags/additional/run_info/pe_runids.txt', 'r') as pe_fh:
    for line in pe_fh:
        pe_runs.add(line.strip())

se_runs = set()
with open('/mfs/gdouglas/projects/ocean_mags/additional/run_info/se_runids.txt', 'r') as se_fh:
    for line in se_fh:
        se_runs.add(line.strip())

with open('/mfs/gdouglas/projects/ocean_mags/additional_round2/samples_formerly_excluded.tsv', 'r') as samples_fh:
    samples_fh.readline()
    for line in samples_fh:
        sample_id, run_ids_raw = line.strip().split('\t')
        run_ids = run_ids_raw.split(',')
        if sample_id in sample_to_runs.keys():
            sys.exit("Error: sample ID present multiple times: " + sample_id)
        else:
            sample_to_runs[sample_id] = run_ids
            sample_to_folder[sample_id] = "/mfs/gdouglas/projects/ocean_mags/additional/filtered"
            
            run_seqtypes = set()
            for run_id in run_ids:
                if run_id in pe_runs:
                    run_seqtypes.add('paired')
                elif run_id in se_runs:
                    run_seqtypes.add('single')
                else:
                    sys.exit("Error: run ID not in PE or SE samples: " + run_id)

            if len(run_seqtypes) > 1:
                sys.exit("Error: multiple seqtypes for sample: " + sample_id)
            
            sample_to_seqtype[sample_id] = list(run_seqtypes)[0]

            # Also add these runs to tmp_runs.
            for run_id in run_ids:
                tmp_runs.add(run_id)

# Then read through additional runs that were totally excluded last time.
missed_runs = set()
with open('/mfs/gdouglas/projects/ocean_mags/additional_round2/additional_run_ids/OceanDNA_nonTara_run_ids.txt', 'r') as missed_fh:
    for line in missed_fh:
        missed_runs.add(line.strip())

# Make sure no runs intersect with the runs that were excluded last time.
if len(tmp_runs.intersection(missed_runs)) > 0:
    sys.exit("Error: runs in both sets.")

with open('/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv', 'r') as ocean_dna_fh:
    col_to_i = {}
    header = ocean_dna_fh.readline().strip().split('\t')
    for i, col in enumerate(header):
        col_to_i[col] = i
    for line in ocean_dna_fh:
        line = line.strip().split('\t')
        sample_id = line[col_to_i['sample_name']]
        run_ids = line[col_to_i['sra_run']].split(',')
        run_ids_set = set(run_ids)

        # Check if any missed_runs in run_ids.
        if len(run_ids_set.intersection(missed_runs)) > 0:
            for run_id in run_ids:
                if run_id not in missed_runs:
                    sys.exit("Error: run_id not in missed set, but corresponds to this sample...") 
        else:
            continue

        if sample_id in sample_to_runs.keys():
            sys.exit("Error: sample ID present multiple times: " + sample_id)
        else:
            sample_to_runs[sample_id] = run_ids
            sample_to_folder[sample_id] = "/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE"
            sample_to_seqtype[sample_id] = 'paired'

num_unique_reps = 0
num_multi_reps = 0
ln_cmds = []
cat_cmds = []

for sample_id in sample_to_runs.keys():
    runs = sample_to_runs[sample_id]
    if len(runs) == 0:
        sys.exit("Error: no runs for " + sample_id + " ?")
    elif len(runs) == 1:
        num_unique_reps += 1
        RUNID = runs[0]

        if sample_to_seqtype[sample_id] == 'paired':
            exp_R1 = sample_to_folder[sample_id] + '/' + RUNID + '_1.fastq.gz'
            exp_R2 = sample_to_folder[sample_id] + '/' + RUNID + '_2.fastq.gz'

            out_R1 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_nonTara/' + sample_id + '_1.fastq.gz'
            out_R2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_nonTara/' + sample_id + '_2.fastq.gz'

            if not os.path.exists(exp_R1) or not os.path.exists(exp_R2):
                sys.exit("Expected files not found: " + exp_R1 + " " + exp_R2)
            else:
                ln_cmds.append("ln -s " + exp_R1 + " " + out_R1)
                ln_cmds.append("ln -s " + exp_R2 + " " + out_R2)
        
        elif sample_to_seqtype[sample_id] == 'single':
            exp_SE = sample_to_folder[sample_id] + '/' + RUNID + '.fastq.gz'

            out_SE = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_SE_final_nonTara/' + sample_id + '.fastq.gz'

            if not os.path.exists(exp_SE):
                sys.exit("Expected file not found: " + exp_SE)
            else:
                ln_cmds.append("ln -s " + exp_SE + " " + out_SE)
            
        else:
            sys.exit("Error: unexpected seqtype: " + sample_to_seqtype[sample_id])

    elif len(runs) > 1:
        num_multi_reps += 1

        cat_cmd_R1 = "cat "
        cat_cmd_R2 = "cat "
        cat_cmd_SE = "cat "

        if sample_to_seqtype[sample_id] == 'paired':

            for RUNID in runs:
                exp_R1 = sample_to_folder[sample_id] + '/' + RUNID + '_1.fastq.gz'
                exp_R2 = sample_to_folder[sample_id] + '/' + RUNID + '_2.fastq.gz'

                if not os.path.exists(exp_R1) or not os.path.exists(exp_R2):
                    sys.exit("Expected files not found: " + exp_R1 + " " + exp_R2)
                else:
                    cat_cmd_R1 += " " + exp_R1
                    cat_cmd_R2 += " " + exp_R2

            out_R1 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_nonTara/' + sample_id + '_1.fastq.gz'
            out_R2 = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_PE_final_nonTara/' + sample_id + '_2.fastq.gz'
            cat_cmd_R1 += " > " + out_R1
            cat_cmd_R2 += " > " + out_R2
            cat_cmds.append(cat_cmd_R1)
            cat_cmds.append(cat_cmd_R2)

        elif sample_to_seqtype[sample_id] == 'single':
            
            for RUNID in runs:
                exp_SE = sample_to_folder[sample_id] + '/' + RUNID + '.fastq.gz'

                if not os.path.exists(exp_SE):
                    sys.exit("Expected file not found: " + exp_SE)
                else:
                    cat_cmd_SE += " " + exp_SE
            
            out_SE = '/mfs/gdouglas/projects/ocean_mags/additional_round2/filtered_SE_final_nonTara/' + sample_id + '.fastq.gz'
            cat_cmd_SE += " > " + out_SE
            cat_cmds.append(cat_cmd_SE)

print("Unique reps:", file=sys.stderr)
print(num_unique_reps, file=sys.stderr)

print("Multi reps:", file=sys.stderr)
print(num_multi_reps, file=sys.stderr)

# Write out commands.
with open("/mfs/gdouglas/projects/ocean_mags/additional_round2/nonTara_ln_cmds.sh", "w") as ln_fh:
    for ln_cmd in ln_cmds:
        print(ln_cmd, file=ln_fh)

with open("/mfs/gdouglas/projects/ocean_mags/additional_round2/nonTara_cat_cmds.sh", "w") as cat_fh:
    for cat_cmd in cat_cmds:
        print(cat_cmd, file=cat_fh)
