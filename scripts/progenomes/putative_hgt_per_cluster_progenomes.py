#!/usr/bin/python3

import argparse
import os
import sys
import gzip
from collections import defaultdict
import hammingdist

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def main():

    parser = argparse.ArgumentParser(

        description='''
Parse aligned FASTAs of sequences from the same cluster.
Identify reciprocal best-hits above 95% (that are at least between genomes of different genera or above).
Also, consider all hits initially, but only output reciprocal best hits that are in focal high-quality genomes.
''',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input',
                        metavar="PATH", type=str,
                        help="Path to folder with paths to input FASTAs.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/clusters/cluster_seqs_aligned/')

    parser.add_argument('-o', '--output',
                        metavar="OUTPUT", type=str,
                        help="Path to output file.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/clusters/putative_hgt_calls.tsv')

    parser.add_argument('-t', '--tax',
                        metavar="TAXTAB", type=str,
                        help="Path to (gzipped) table with taxonomic information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/progenomes_ncbi_taxonomy.tsv.gz')

    parser.add_argument('--highqual_genomes',
                        metavar="IDS", type=str,
                        help="Path to file with high quality genomes (with metagenomics assemblies and optionally, ambiguous genera excluded).",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/clearly_non_mags_genome_ids.txt')
 
    parser.add_argument('-g', '--gene_to_contig',
                        metavar="PATH", type=str,
                        help="Path to (gzipped) table with mapping of gene to genome and contig.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/gene_to_contig.tsv.gz')

    parser.add_argument('-c', '--contig_lengths',
                        metavar="PATH", type=str,
                        help="Path to (gzipped) table with contig lengths (and categorized as long or short, along with whether it's a plasmid).",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/contig_length_breakdown.tsv.gz')

    args = parser.parse_args()

    highqual_genomes = set()
    with open(args.highqual_genomes, 'r') as genomes_fh:
        for line in genomes_fh:
            highqual_genomes.add(line.strip())

    # Read in genome taxonomy.
    tax = defaultdict(dict)
    tax_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']
    tax_col_to_i = {}
    with gzip.open(args.tax, "rb") as tax_fh:
        tax_header = tax_fh.readline().decode("utf-8").strip().split("\t")
        for i, col in enumerate(tax_header):
            tax_col_to_i[col] = i

        for tax_line in tax_fh:
            tax_line = tax_line.decode("utf-8").strip().split("\t")
            taxid = tax_line[tax_col_to_i["taxid"]]
            for tax_level in tax_levels:
                tax[taxid][tax_level] = tax_line[tax_col_to_i[tax_level]]

    gene_to_contig = {}
    with gzip.open(args.gene_to_contig, "rb") as gene_to_contig_fh:
        gene_to_contig_fh.readline()
        for gene_to_contig_line in gene_to_contig_fh:
            gene_to_contig_line = gene_to_contig_line.decode("utf-8").strip().split("\t")
            gene_to_contig[gene_to_contig_line[0]] = gene_to_contig_line[3]

    passing_contigs = set()
    with gzip.open(args.contig_lengths, "rb") as contig_lengths_fh:
        contig_lengths_fh.readline()
        for contig_lengths_line in contig_lengths_fh:
            contig_lengths_line = contig_lengths_line.decode("utf-8").strip().split("\t")
            if contig_lengths_line[2] == 'long' or contig_lengths_line[4] == 'Plasmid':
                passing_contigs.add(contig_lengths_line[1])    

    with open(args.output, 'w') as out_fh:
        print('\t'.join(['gene1', 'gene2', 'gene1_genome', 'gene2_genome', 'highest_tax_diff', 'identity']), file=out_fh)
        in_fastas = os.listdir(args.input)
        for fasta in in_fastas:
            filepath = os.path.join(args.input, fasta)
            seq_ids = []
            seqs = {}
            name = None
            with open(filepath, 'r') as fasta_in:
                for line in fasta_in:
                    line = line.rstrip()
                    if len(line) == 0:
                        continue
                    # If header-line then split by whitespace, take the first element,
                    # and define the sequence name as everything after the ">".
                    if line[0] == ">":
                        name = line[1:]
                        if name in seqs:
                            sys.stderr("Stopping due to duplicated id in file: " + name)
                        seqs[name] = ''
                        seq_ids.append(name)
                    else:
                        seqs[name] += line

            # Make sure all sequences are same length.
            seq_len = len(seqs[list(seqs.keys())[0]])
            for seq in seqs.values():
                if len(seq) != seq_len:
                    sys.exit("Error: Not all sequences are the same length.")

            # Skip if seq_len <= 500.
            if seq_len <= 500:
                continue

            # Get pairwise distances of aligned sequences.
            dist = hammingdist.from_fasta(filepath)

            top_matches = {}
            top_identities = {}
            ambig_seqs = set()
            for i in range(len(seq_ids) - 1):
                for j in range(i + 1, len(seq_ids)):
                    i_id = seq_ids[i]
                    j_id = seq_ids[j]
                    iden = (1 - (dist[i, j] / seq_len)) * 100.0
                    if iden >= 95.0:
                        if i_id not in top_matches or iden > top_identities[i_id]:
                            top_matches[i_id] = j_id
                            top_identities[i_id] = iden
                        elif iden == top_identities[i_id]:
                            ambig_seqs.add(i_id)

                        if j_id not in top_matches or iden > top_identities[j_id]:
                            top_matches[j_id] = i_id
                            top_identities[j_id] = iden
                        elif iden == top_identities[j_id]:
                            ambig_seqs.add(j_id)

            parsed_seqs = set()
            for seq_id in seq_ids:
                if seq_id not in top_matches or seq_id in ambig_seqs or seq_id in parsed_seqs:
                    continue

                best_match = top_matches[seq_id]
                reverse_best_match = top_matches[best_match]

                if seq_id == reverse_best_match and best_match not in ambig_seqs and best_match not in parsed_seqs:
                    parsed_seqs.add(best_match)
                    parsed_seqs.add(seq_id)

                    seq_id_split = seq_id.split(".")
                    gene1_genome = 'g' + seq_id_split[0] + "_" + seq_id_split[1]
                    gene1_taxid = seq_id_split[0]

                    best_match_split = best_match.split(".")
                    gene2_genome = 'g' + best_match_split[0] + "_" + best_match_split[1]
                    gene2_taxid = best_match_split[0]

                    if gene1_genome == gene2_genome or gene1_taxid == gene2_taxid:
                        continue

                    if gene1_genome not in highqual_genomes or gene2_genome not in highqual_genomes:
                        continue

                    gene1_contig = gene_to_contig[seq_id]
                    gene2_contig = gene_to_contig[best_match]
                    if gene1_contig not in passing_contigs or gene2_contig not in passing_contigs:
                        continue

                    diff_found = False
                    for tax_level in tax_levels:
                        if tax[gene1_taxid][tax_level] != tax[gene2_taxid][tax_level]:
                            highest_tax_diff = tax_level
                            diff_found = True
                            break

                    if not diff_found:
                        print("Warning: No taxonomic difference found between genomes " + gene1_genome + " and " + gene2_genome + ".", file=sys.stderr)
                        continue

                    if highest_tax_diff in ['Species', 'Strain']:
                        continue

                    print("\t".join([seq_id, best_match, gene1_genome, gene2_genome, highest_tax_diff, str(top_identities[seq_id])]), file=out_fh)


if __name__ == '__main__':
    main()
