#!/usr/bin/python3

import os
import gzip

# Get mapfile of genes to Panaroo gene families.

panaroo_out_dir = "/mfs/gdouglas/projects/ocean_mags/species_DTL_analyses/panaroo_out"
panaroo_out_folders = os.listdir(panaroo_out_dir)

print("gene_id\tpanaroo_sp_and_gene_family")

for panaroo_subfolder in panaroo_out_folders:
    gene_data_file = os.path.join(panaroo_out_dir, panaroo_subfolder, "gene_presence_absence.csv.gz")

    with gzip.open(gene_data_file, 'rt') as gene_data_fh:
        headerline = gene_data_fh.readline().strip().split(",")

        for line in gene_data_fh:
            line = line.strip().split(",")
            gene_family = line[0]

            for i in range(3, len(line)):
                gene_id_raw = line[i]
                if gene_id_raw == '':
                    continue
                
                if ';' in gene_id_raw:
                    gene_ids = gene_id_raw.split(";")
                else:
                    gene_ids = [gene_id_raw]

                for gene_id in gene_ids:
                    print(f"{gene_id}\t{panaroo_subfolder},{gene_family}")
