#!/usr/bin/python3

import os
import sys
from collections import defaultdict

# Parse folder with FASTP-filtered FASTQs and:
# (1) Keep only one sample (with the highest % of prokaryotic DNA) for cases where samples share all the same collection information (but just differing filters!).
# (2) For run IDs that are the sole runs for a sample, produce symolic links to renamed FASTQs that include the sample ID instead.
# (3) For run IDs that are one of multiple of the same sample, produce cat commands to concatenate the FASTQ(s) for that sample.

# First, read through metadata table to figure out which samples are redundant based on metadata.
redundant_samples = defaultdict(list)
with open("/mfs/gdouglas/projects/water_mags/additional/OceanDNA_supp_metadata/subset_tab.tsv", "r") as metadata_fh:
    metadata_col = metadata_fh.readline().strip().split("\t")
    metadata_col_to_i = {col: i for i, col in enumerate(metadata_col)}

    for line in metadata_fh:
        line = line.strip().split("\t")
        sample_id = line[metadata_col_to_i["sample_name"]]
        collection_date = line[metadata_col_to_i["collection_date"]]
        depth = line[metadata_col_to_i["depth"]]
        latitude = line[metadata_col_to_i["latitude"]]
        # Note that the typo "longigute" was in the original metadata table.
        longitude = line[metadata_col_to_i["longigute"]]
        percent_prok = line[metadata_col_to_i["kaiju_percent_of_Bacteria"]]

        collection_info = (collection_date, depth, latitude, longitude)
        sample_id_info = (sample_id, float(percent_prok))
        redundant_samples[collection_info].append(sample_id_info)

nonindependent_samples_to_exclude = set()
for collection_info, sample_info_list in redundant_samples.items():
    if len(sample_info_list) > 1:
        best_sample = sample_info_list[0][0]
        highest_prok = sample_info_list[0][1]
        all_samples = set([best_sample])
        for i in range(1, len(sample_info_list)):
            all_samples.add(sample_info_list[i][0])
            if sample_info_list[i][1] > highest_prok:
                best_sample = sample_info_list[i][0]
                highest_prok = sample_info_list[i][1]

        nonindependent_samples_to_exclude.update(all_samples - set([best_sample]))

# Then read through all filtered FASTQ names to get map of run ID to FASTQ files.
additional_ids = set()
runid_to_files = defaultdict(set)
extra_file_list = os.listdir("/mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered")
for file in extra_file_list:
    runid = file.split("_")[0]
    runid = runid.replace(".fastq.gz", "")
    additional_ids.add(runid)
    runid_to_files[runid].add("/mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered/" + file)

file_list = os.listdir("/mfs/gdouglas/projects/water_mags/additional/filtered")
for file in file_list:
    runid = file.split("_")[0]
    runid = runid.replace(".fastq.gz", "")
    if runid in additional_ids:
        continue
    runid_to_files[runid].add("/mfs/gdouglas/projects/water_mags/additional/filtered/" + file)

parsed_runs = set()

ln_cmds = []
cat_cmds = []

se_samples = []
pe_samples = []
with open("/mfs/gdouglas/projects/water_mags/additional/OceanDNA_supp_metadata/subset_tab.tsv", "r") as metadata_fh:
    metadata_fh.readline()

    for line in metadata_fh:
        line = line.strip().split("\t")
        sample_id = line[metadata_col_to_i["sample_name"]]

        if sample_id in nonindependent_samples_to_exclude:
            continue

        sample_runids = line[metadata_col_to_i["sra_run"]].split(',')
        sample_runids_set = set(sample_runids)
        # Skip if no additional ids intersect with these run IDs.
        if len(sample_runids_set & additional_ids) == 0:
            continue

        # Create symbolic link for cases where a single run ID corresponds to sample.
        if len(sample_runids) == 1:
            runid = sample_runids[0]
            if runid in parsed_runs:
                sys.exit("Error: run ID already parsed: " + runid)
            parsed_runs.add(runid)

            # Write out symbolic link command.
            run_fastqs = runid_to_files[runid]

            run_fastqs_base = []
            for run_fastq in run_fastqs:
                run_fastq_base = os.path.basename(run_fastq)
                run_fastqs_base.append(run_fastq_base)
            
            if len(run_fastqs) == 1:
                exp_se = runid + ".fastq.gz"
                if exp_se != run_fastqs_base[0]:
                    sys.exit("Error: unexpected filename for: " + runid)
                exp_se = list(run_fastqs)[0]
                outfile_se = sample_id + ".fastq.gz"
                ln_cmds.append("ln -s " + exp_se + " /mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean/" + outfile_se)
                se_samples.append(sample_id)
            elif len(run_fastqs) == 2:
                exp_r1 = runid + "_1.fastq.gz"
                exp_r2 = runid + "_2.fastq.gz"
                if exp_r1 not in run_fastqs_base or exp_r2 not in run_fastqs_base:
                    sys.exit("Error: unexpected filenames for: " + runid)

                run_fastqs = list(run_fastqs)
                if exp_r1 in run_fastqs[0]:
                    exp_r1 = run_fastqs[0]
                    exp_r2 = run_fastqs[1]
                else:
                    exp_r1 = run_fastqs[1]
                    exp_r2 = run_fastqs[0]

                outfile_r1 = sample_id + "_1.fastq.gz"
                outfile_r2 = sample_id + "_2.fastq.gz"
                ln_cmds.append("ln -s " + exp_r1 + " /mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean/" + outfile_r1)
                ln_cmds.append("ln -s " + exp_r2 + " /mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean/" + outfile_r2)
                pe_samples.append(sample_id)
            else:
                sys.exit("Error: unexpected number of files for run ID: " + runid)
        elif len(sample_runids) > 1:
            se_fastqs = []
            r1_fastqs = []
            r2_fastqs = []

            # Write out cat command for cases where multiple run IDs correspond to sample.
            for runid in sample_runids:
                if runid in parsed_runs:
                    sys.exit("Error: run ID already parsed: " + runid)
                parsed_runs.add(runid)

                run_fastqs = runid_to_files[runid]

                run_fastqs_base = []
                for run_fastq in run_fastqs:
                    run_fastq_base = os.path.basename(run_fastq)
                    run_fastqs_base.append(run_fastq_base)

                if len(run_fastqs) == 1:
                    exp_se = runid + ".fastq.gz"
                    if exp_se != list(run_fastqs_base)[0]:
                        sys.exit("Error: unexpected filename for: " + runid)
                    exp_se = list(run_fastqs)[0]
                    se_fastqs.append(exp_se)
                elif len(run_fastqs) == 2:
                    exp_r1 = runid + "_1.fastq.gz"
                    exp_r2 = runid + "_2.fastq.gz"
                    if exp_r1 not in run_fastqs_base or exp_r2 not in run_fastqs_base:
                        print(exp_r1, exp_r2, run_fastqs_base)
                        sys.exit("Error: unexpected filenames for: " + runid)
                    run_fastqs = list(run_fastqs)
                    if exp_r1 in run_fastqs[0]:
                        exp_r1 = run_fastqs[0]
                        exp_r2 = run_fastqs[1]
                    else:
                        exp_r1 = run_fastqs[1]
                        exp_r2 = run_fastqs[0]
                    r1_fastqs.append(exp_r1)
                    r2_fastqs.append(exp_r2)
                else:
                    sys.exit("Error: unexpected number of files for run ID: " + runid + ": " + str(len(run_fastqs)))
            
            if len(se_fastqs) > 0 and (len(r1_fastqs) > 0 or len(r2_fastqs) > 0):
                sys.exit("Error: unexpected mix of SE and PE files for sample: " + sample_id)

            elif len(se_fastqs) > 0:
                cat_cmds.append("cat " + " ".join(se_fastqs) + " > /mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean/" + sample_id + ".fastq.gz")
                se_samples.append(sample_id)
            
            elif len(r1_fastqs) > 0:
                if len(r1_fastqs) != len(r2_fastqs):
                    sys.exit("Error: unexpected number of R1 and R2 files for sample: " + sample_id)
                cat_cmds.append("cat " + " ".join(r1_fastqs) + " > /mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean/" + sample_id + "_1.fastq.gz")
                cat_cmds.append("cat " + " ".join(r2_fastqs) + " > /mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean/" + sample_id + "_2.fastq.gz")
                pe_samples.append(sample_id)
            else:
                sys.exit("No FASTQs for sample??: " + sample_id)

# Write out symbolic link and cat commands.
with open("/mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean_ln_cmds.sh", "w") as ln_fh:
   
   for ln_cmd in ln_cmds:
       ln_fh.write(ln_cmd + "\n")

with open("/mfs/gdouglas/projects/water_mags/additional/failed_redownloaded/filtered_clean_cat_cmds.sh", "w") as cat_fh:
   for cat_cmd in cat_cmds:
       cat_fh.write(cat_cmd + "\n")

# # Write out PE and SE sample IDs.
# with open("/mfs/gdouglas/projects/water_mags/additional/run_info/se_sampleids.txt", "w") as se_fh:
#     for sample_id in sorted(list(se_samples)):
#         se_fh.write(sample_id + "\n")

# with open("/mfs/gdouglas/projects/water_mags/additional/run_info/pe_sampleids.txt", "w") as pe_fh:
#     for sample_id in sorted(list(pe_samples)):
#         pe_fh.write(sample_id + "\n")
