#!/usr/bin/python3

import os
import sys
from collections import defaultdict

# Quick script to figure out which run IDs correspond to SE vs PE data.

file_list = os.listdir("/mfs/gdouglas/projects/water_mags/additional/raw")

runid_to_files = defaultdict(set)

for file in file_list:
    runid = file.split("_")[0]
    runid = runid.replace(".fastq.gz", "")

    runid_to_files[runid].add(file)

se = set()
pe = set()
extra_files = set()

for runid, files in runid_to_files.items():
    if len(files) == 3:
        # There are some cases where there is an additional file with se or technical reads.
        # These seem to be much smaller files, so will print these out for further investigation
        # (and likely deletion from the raw directory), and treat the run ID as PE.
        r1_expected = runid + "_1.fastq.gz"
        r2_expected = runid + "_2.fastq.gz"
        r3_expected = runid + ".fastq.gz"
        if r1_expected in files and r2_expected in files and r3_expected in files:
            pe.add(runid)
            extra_files.add(r3_expected)
        else:
            sys.exit("Error: unexpected files for run ID: " + runid)

    elif len(files) == 1:
        se.add(runid)
    elif len(files) == 2:
        pe.add(runid)
    else:
        sys.exit("Error: unexpected number of files for run ID: " + runid)

# Write out run IDs and extra filenames.
with open("/mfs/gdouglas/projects/water_mags/additional/run_info/se_runids.txt", "w") as f:
    for runid in sorted(list(se)):
        f.write(runid + "\n")

with open("/mfs/gdouglas/projects/water_mags/additional/run_info/pe_runids.txt", "w") as f:
    for runid in sorted(list(pe)):
        f.write(runid + "\n")

with open("/mfs/gdouglas/projects/water_mags/additional/run_info/extra_files.txt", "w") as f:
    for filename in extra_files:
        f.write(filename + "\n")
