#!/usr/bin/python3

import gzip
import os
import sys
from itertools import combinations

# Get pairwise HGT inferences for every pair of genomes in the summarized output.
def main():

    # Get mapping of taxa to genome names.
    with gzip.open("/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_to_taxon.tsv.gz", 'rt') as taxon_map_fh:
        taxon_map = dict()
        for line in taxon_map_fh:
            line_split = line.split()
            taxon_map[line_split[0]] = line_split[1]

    print("taxa_combo\tspecies\ttaxon1\ttaxon2\thgt_count")

    for species in sorted(os.listdir("/mfs/gdouglas/projects/ocean_mags/water_mag_analysis/species_DTL_analyses/homer_rangerdtl_summaries")):
        
        # Read in all genomes used for this species.
        genome_map_file = "/mfs/gdouglas/projects/ocean_mags/water_mag_analysis/species_DTL_analyses/homer_prep_map_only/" + species + "/map_out/genome_ids.tsv"
        genomes = []
        with open(genome_map_file, 'r') as genome_map_fh:
            for genome_line in genome_map_fh:
                genome_line = genome_line.strip().split()
                genomes.append(taxon_map[genome_line[0]])

        # Then initialize the HGT count for every taxa combo to 0.
        hgt_count = {}
        genomes = sorted(genomes)
        taxa_combinations = list(combinations(genomes, 2))
        for taxa_combo in taxa_combinations:
            hgt_count[','.join(taxa_combo)] = 0

        with open("/mfs/gdouglas/projects/ocean_mags/water_mag_analysis/species_DTL_analyses/homer_rangerdtl_summaries/" + species + "/transfers.tsv", 'r') as transfers_fh:
            colmap = dict()
            col_names_list = transfers_fh.readline().split()
            for i in range(len(col_names_list)):
                colmap[col_names_list[i]] = i
            
            for line in transfers_fh:
                line_split = line.split()

                donor = line_split[colmap["most.freq.donor"]]
                recipient = line_split[colmap["most.freq.recipient"]]
                
                if donor not in taxon_map:
                    if len(donor) < 5:
                        continue
                    else:
                        sys.exit("Error: " + donor + " not found in taxon map.")
                elif recipient not in taxon_map:
                    if len(recipient) < 5:
                        continue
                    else:
                        sys.exit("Error: " + recipient + " not found in taxon map.")

                if line_split[colmap["most.freq.donor.instances"]] == "100" and line_split[colmap["most.freq.recipient.instances"]] == "100":
                    donor_taxon = taxon_map[donor]
                    recipient_taxon = taxon_map[recipient]
                    taxa_combo = ','.join(sorted([donor_taxon, recipient_taxon]))
                    hgt_count[taxa_combo] += 1

        for taxa_combo in hgt_count:
            taxa_split = taxa_combo.split(',')
            print(taxa_combo + "\t" + species + "\t" + taxa_split[0] + "\t" + taxa_split[1] + "\t" + str(hgt_count[taxa_combo]))


if __name__ == '__main__':
    main()
