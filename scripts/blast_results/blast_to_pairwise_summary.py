#!/usr/bin/python3

from collections import defaultdict
import gzip
import os
import sys

# Read in genome taxonomy.
tax = defaultdict(dict)
tax_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]

# Read in mappings of genome IDs to taxonomy.
tax_col_to_i = {}
with gzip.open("/mfs/gdouglas/projects/water_mags/mapfiles/MAG_taxa_breakdown.tsv.gz", "rb") as tax_fh:
    tax_header = tax_fh.readline().decode("utf-8").strip().split("\t")
    for i, col in enumerate(tax_header):
        tax_col_to_i[col] = i

    for tax_line in tax_fh:
        tax_line = tax_line.decode("utf-8").strip().split("\t")
        taxon_id = tax_line[tax_col_to_i["Taxon_ID"]]
        for tax_level in tax_levels:
            tax[taxon_id][tax_level] = tax_line[tax_col_to_i[tax_level]]

# Read in all BLAST gene count bedfiles.
gene_counts = {}
count_bedfiles = os.listdir("/mfs/gdouglas/projects/water_mags/blast_output/intersections/blast_hit_gene_counts")
for bedfile in count_bedfiles:
    bedfile_id = bedfile.replace("Hits_", "").replace("_genes.bed", "")
    gene_counts[bedfile_id] = {}
    file_percent = int(bedfile_id.split("_")[-1])
    with open("/mfs/gdouglas/projects/water_mags/blast_output/intersections/blast_hit_gene_counts/" + bedfile, "r") as bed_fh:
        for bed_line in bed_fh:
            bed_line = bed_line.strip().split("\t")
            gene_counts[bedfile_id][bed_line[3]] = int(bed_line[6])
            line_percent = float(bed_line[4])
            if file_percent == 95:
                if line_percent < 95 or line_percent >= 99:
                    sys.exit("Error: Hit percent is less than 95 or higher than 99.")
            elif file_percent == 99:
                if line_percent < 99:
                    sys.exit("Error: Hit percent is less than 99.")

# List all BLAST filtered files.
blast_files = os.listdir("/mfs/gdouglas/projects/water_mags/blast_output/dedup_blast_out")
print("Taxa_combo\tHighest_tax_diff\t95_hit_count\t99_hit_count\tboth_hit_count\t95_gene_count\t99_gene_count\tboth_gene_count")
for blast_file in blast_files:
    hit_counts_95 = 0
    gene_counts_95 = 0
    hit_counts_99 = 0
    gene_counts_99 = 0
    hit_counts_both = 0
    gene_counts_both = 0
    highest_tax_diff = "None found"

    blast_file_split = blast_file.split(".")
    taxon1 = blast_file_split[1]
    taxon2 = blast_file_split[3]
    taxa_combo = ",".join(sorted([taxon1, taxon2]))

    diff_found = False
    for tax_level in tax_levels:
        if tax[taxon1][tax_level] != tax[taxon2][tax_level]:
            highest_tax_diff = tax_level
            diff_found = True
            break
    if not diff_found:
        sys.exit("Error: No taxonomic difference found between genomes " + taxon1 + " and " + taxon2 + ".")

    with gzip.open("/mfs/gdouglas/projects/water_mags/blast_output/dedup_blast_out/" + blast_file, "rb") as blast_fh:
        for blast_line in blast_fh:
            blast_line = blast_line.decode("utf-8").strip().split("\t")
            identity = float(blast_line[2])
            if len(blast_line) != 16:
                sys.exit("Error: BLAST output does not have 16 columns.")

            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
            scaffold1 = blast_line[0]
            scaffold2 = blast_line[1]
            qstart = blast_line[6]
            qend = blast_line[7]
            sstart = blast_line[8]
            send = blast_line[9]

            query_id = "query_" + scaffold1 + "||" + qstart + "||" + qend + "||" + scaffold2 + "||" + sstart + "||" + send
            subject_id = "subject_" + scaffold1 + "||" + qstart + "||" + qend + "||" + scaffold2 + "||" + sstart + "||" + send

            if identity >= 99:
                hit_counts_99 += 1
                hit_compare1 = gene_counts[highest_tax_diff + "_99"][query_id]
                hit_compare2 = gene_counts[highest_tax_diff + "_99"][subject_id]
                gene_counts_99 += (hit_compare1 + hit_compare2) / 2
                gene_counts_both += (hit_compare1 + hit_compare2) / 2
            elif identity >= 95 and identity < 99:
                hit_counts_95 += 1
                hit_compare1 = gene_counts[highest_tax_diff + "_95"][query_id]
                hit_compare2 = gene_counts[highest_tax_diff + "_95"][subject_id]
                gene_counts_95 += (hit_compare1 + hit_compare2) / 2
                gene_counts_both += (hit_compare1 + hit_compare2) / 2
            else:
                sys.exit("Error: Identity outside of expected range.")

            hit_counts_both += 1

    if hit_counts_both > 0:
        print(taxa_combo + "\t" + highest_tax_diff + "\t" + str(hit_counts_95) + "\t" + str(hit_counts_99) + "\t" + str(hit_counts_both) + "\t" + str(gene_counts_95) + "\t" + str(gene_counts_99) + "\t" + str(gene_counts_both))

print("Done.", file=sys.stderr)
