#!/usr/bin/python3

import os
import sys
from collections import defaultdict

# Parse folders with FASTP-filtered FASTQs and:
# (1) For run IDs that are the sole runs for a sample, produce symolic links to renamed FASTQs that include the sample ID instead.
# (2) For run IDs that are one of multiple of the same sample, produce cat commands to concatenate the FASTQ(s) for that sample.

# Read through Tara ocean metadata tables first.
# This RNA-seq table needs to be run separately: 'PRJEB6608'
TARA_bioprojects = ['PRJEB1787', 'PRJEB9740']
TARA_sample_to_run = defaultdict(list)
parsed_runs = set()
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
            sample_ID = line[sample_col_i]

            if run_ID in parsed_runs:
                sys.exit("Error: run ID already parsed: " + run_ID)
            else:
                parsed_runs.add(run_ID)
            TARA_sample_to_run[sample_ID].append(run_ID)

Tara_samples_from_csv = set(list(TARA_sample_to_run.keys()))
Tara_samples_in_Table = set()

if 'ERS477931' in Tara_samples_from_csv:
    print("Yes!")


with open("/mfs/gdouglas/projects/ocean_mags/metadata/OceanDNA_supp_metadata/Supp_File_S1_water_samples.tsv", "r") as Tab_fh:
    metadata_col = Tab_fh.readline().strip().split("\t")
    metadata_col_to_i = {col: i for i, col in enumerate(metadata_col)}

    for line in Tab_fh:
        line = line.strip().split("\t")
        sample_id = line[metadata_col_to_i["sample_name"]]

        if sample_id not in Tara_samples_from_csv:
            continue

        Tara_samples_in_Table.add(sample_id)

        sra_runs = sorted(line[metadata_col_to_i["sra_run"]].split(','))

        if sra_runs != sorted(TARA_sample_to_run[sample_id]):

            if len(set(sra_runs) - set(sorted(TARA_sample_to_run[sample_id]))) > 0:
                sys.exit("Error: SRA run ID mismatch for sample ID: " + sample_id)

if Tara_samples_in_Table != Tara_samples_from_csv:
    print(len(Tara_samples_in_Table))
    print(len(Tara_samples_from_csv))
    print(len(Tara_samples_in_Table - Tara_samples_from_csv))
    print(len(Tara_samples_from_csv - Tara_samples_in_Table))
    sys.exit("Error: sample ID mismatch between metadata tables and Tara ocean metadata tables.")
