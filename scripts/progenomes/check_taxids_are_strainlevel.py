#!/usr/bin/env python3
import sys
from Bio import Entrez

###from collections import defaultdict
##taxonomy = {}
##levels_of_interest = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']

# Quick check that taxids in proGenomes labels are all at the strain level (so that it makes sense to group scaffolds based on these IDs).
Entrez.email = "gavinmdouglas@gmail.com"

with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/all_taxids.txt') as taxid_fh:
    taxids = [line.strip() for line in taxid_fh if line.strip()]

# Process in batches to be efficient with API calls
batch_size = 100
non_strain_taxids = []

for i in range(0, len(taxids), batch_size):
    batch = taxids[i:i + batch_size]

    # Fetch taxonomy information
    handle = Entrez.efetch(db="taxonomy", id=batch)
    records = Entrez.read(handle)

    # Also track IDs that could not be found.
    found_taxids = [record['TaxId'] for record in records]
    missing_taxids = set(batch) - set(found_taxids)
    if missing_taxids:
        print(f"Could not find taxonomy information for the following taxids: {missing_taxids}")

    for record in records:
        taxid = record['TaxId']
        taxid_str = 'taxid' + str(taxid)
        rank = record['Rank']
        name = record['ScientificName']

        if rank != 'strain':
            non_strain_taxids.append((taxid, rank, name))

# Report results
if non_strain_taxids:
    print("Found taxids that are not strain level:")
    for taxid, rank, name in non_strain_taxids:
        print(f"Taxid: {taxid}, Rank: {rank}, Name: {name}")
    sys.exit(1)
else:
    print("All taxids are strain level classifications.")
    sys.exit(0)

