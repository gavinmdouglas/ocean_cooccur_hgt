#!/usr/bin/python3

import os
import sys

# Quick check that all downloaded IDs correspond to those in globus batch file.

parsed_ids = set()
with open("/mfs/gdouglas/projects/water_mags/additional/globus_batch.txt", "r") as batch_fh:
    for line in batch_fh:
        runid = line.split()[2].split("/")[-1]
        if runid in parsed_ids:
            sys.exit("Already parsed this ID: %s" % runid)
        else:
            parsed_ids.add(runid)

obs_ids = set()
with open("/mfs/gdouglas/projects/water_mags/additional/run_ids.txt", "r") as obs_ids_fh:
    for line in obs_ids_fh:
        if line.strip() not in parsed_ids:
            print("Not in batch file: %s" % line.strip())
        else:
            obs_ids.add(line.strip())

if obs_ids != parsed_ids:
    sys.exit("Mismatch between batch and obs IDs")
else:
    print("All IDs match")
