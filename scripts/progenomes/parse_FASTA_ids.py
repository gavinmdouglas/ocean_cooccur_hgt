#!/usr/bin/python3

import gzip
from collections import defaultdict

genomes_to_contigs = defaultdict(list)

all_biosample_ids = set()

all_taxids = set()

with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/representatives.aquatic.contigs.headers.txt.gz', 'rt') as fasta_header_fh:
    for line in fasta_header_fh:
        if line.startswith('>'):
            line = line[1:].strip()
            contig_subset = line.split(' ')[0]
            contig_subset_split = contig_subset.split('.')
            genome_id = contig_subset_split[0] + '.' + contig_subset_split[1]
            all_biosample_ids.add(contig_subset_split[1])
            all_taxids.add(contig_subset_split[0])
            genomes_to_contigs[genome_id].append(line)

print('Number of genomes:')
print(len(genomes_to_contigs.keys()))

num_less_than_x = 0
num_x_or_more = 0
for genome_id, contigs in genomes_to_contigs.items():
    if len(contigs) < 10:
        num_less_than_x += 1
        if 'MAG' in genome_id:
            print('Seems complete, but has \"MAG\" in the ID:')
            print(genome_id)
    else:
        num_x_or_more += 1

print('Number of genomes with less than 10 contigs:')
print(num_less_than_x)

print('Number of genomes with 10 or more contigs:')
print(num_x_or_more)

def genomes_w_strings(strings_to_search: list):
    genome_to_match = {}
    num_w_strings = 0
    num_wout_strings = 0
    for genome_id, contigs in genomes_to_contigs.items():
        has_string = False
        for contig in contigs:
            for string in strings_to_search:
                if string in contig:
                    has_string = True
                    break
        if has_string:
            genome_to_match[genome_id] = 'Yes'
            num_w_strings += 1
        else:
            genome_to_match[genome_id] = 'No'
            num_wout_strings += 1
    return genome_to_match, num_w_strings, num_wout_strings

genome_to_MAG, num_w_MAG, num_wout_MAG = genomes_w_strings(['MAG'])
genome_to_mag, num_w_mag, num_wout_mag = genomes_w_strings(['mag'])
genome_to_shotgun, num_w_shotgun, num_wout_shotgun = genomes_w_strings(['shotgun', 'Shotgun', 'SHOTGUN'])

print('Number of genomes with \"MAG\" in the ID:')
print(num_w_MAG)

print('Number of genomes without \"MAG\" in the ID:')
print(num_wout_MAG)


print('Number of genomes with \"mag\" in the ID:')
print(num_w_mag)

print('Number of genomes without \"mag\" in the ID:')
print(num_wout_mag)

print('Number of genomes with \"shotgun\" in the ID:')
print(num_w_shotgun)

print('Number of genomes without \"shotgun\" in the ID:')
print(num_wout_shotgun)
    
    


with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_or_chrom_per_genome.txt', 'w') as out_fh:
    out_fh.write('genome\tcontig_count\tMAG_in_name\tmag_in_name\tshotgun_in_name\n')
    for genome_id, contigs in genomes_to_contigs.items():
        out_fh.write(genome_id + '\t' + str(len(contigs)) + '\t' + genome_to_MAG[genome_id] + '\t' + genome_to_mag[genome_id] + '\t' + genome_to_shotgun[genome_id] + '\n')

with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/all_biosample_ids.txt', 'w') as out_fh:
    for biosample_id in sorted(list(all_biosample_ids)):
        out_fh.write(biosample_id + '\n')

with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/all_taxids.txt', 'w') as out_fh:
    for taxid in sorted(list(all_taxids)):
        out_fh.write(taxid + '\n')
