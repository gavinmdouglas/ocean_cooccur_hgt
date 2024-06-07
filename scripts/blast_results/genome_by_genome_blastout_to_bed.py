#!/usr/bin/python3

import argparse
from os import listdir
from os.path import isfile, join
import sys
import pandas as pd
import gzip

def main():

    parser = argparse.ArgumentParser(

        description='Read through directory full of gzipped BLAST output tables.  '
                    'Importantly, this script assumes that the BLAST output is in this order: '
                    'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore, '
                    'and these columns can optionally be at the end: qlen and slen. '
                    'Write out different bedfile for subset of hits that are between different strains, '
                    'different species, different genera, etc. \n'
                    'Will output separate set of bedfiles for hits in the range of >= 95 and < 99 identity as well as >=99% percent identity.',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--blast_results", metavar="FRAGMENT_BED", type=str,
                        help="Folder containing all BLAST results ", required=True)

    parser.add_argument("--out_prefix", metavar="OUTPUT", type=str,
                        help="Output prefix", required=True)

    parser.add_argument("--taxonomy", metavar="TAXA", type=str,
                        help="Taxonomy table per genome", required=True)

    parser.add_argument("--id_map", metavar="Genome2Taxon", type=str,
                        help="Mapfile mapping raw genome names to Taxon numbers", required=True)

    parser.add_argument("--contig2genome", metavar="Contig2genome", type=str,
                        help="Mapfile of scaffold/contig to genomes", required=True)

    args = parser.parse_args()

    # Cleaned-up taxonomic levels per genome.
    genome_taxonomy = pd.read_csv(args.taxonomy,
                                  index_col = 0,
                                  header = 0,
                                  sep = '\t')

    # Map of contig ids to genome ids (e.g., BATS_SAMN07137064_METAG-scaffold_10025 BATS_SAMN07137064_METAG_ADEDPHOH)
    contig2genome = pd.read_csv(args.contig2genome,
                                index_col = 0,
                                header = None,
                                sep = ' ',
                                names = ['contig', 'genome'])

    # Map of original genome ids to simplified taxon ids (e.g., BATS_SAMN07137064_METAG_BDCGNAHA        Taxa_1)
    orig2simple = pd.read_csv(args.id_map,
                              index_col = 0,
                              header = None,
                              sep = '\t',
                              names = ['original', 'taxon'])

    blast_results = [join(args.blast_results, f) for f in listdir(args.blast_results) if isfile(join(args.blast_results, f))]

    identity_cutoffs = ['95', '99']

    tax_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']

    for identity_cutoff in identity_cutoffs:

        out_filehandles = dict()
        for tax_level in tax_levels:
            out_filehandles[tax_level] = open(args.out_prefix + '_' + tax_level + '_' + identity_cutoff + '.bed', 'w')

        for blast_result in blast_results:

            # Read through BLAST files and print out lines where the matched contigs differ by the relevant taxonomic label.
            with gzip.open(blast_result, mode="rt") as in_blast_filehandle:

                for line in in_blast_filehandle:

                    line_split = line.split('\t')

                    query_contig = line_split[0]
                    query_start = int(line_split[6])
                    query_end = int(line_split[7])

                    subject_contig = line_split[1]
                    subject_start = int(line_split[8])
                    subject_end = int(line_split[9])

                    query_genome = contig2genome.loc[query_contig, 'genome']
                    subject_genome = contig2genome.loc[subject_contig, 'genome']

                    percent_identity = float(line_split[2])

                    if identity_cutoff == '95' and (percent_identity >= 99.0 or percent_identity < 95.0):
                        continue
                    elif identity_cutoff == '99' and percent_identity < 99.0:
                        continue

                    if query_genome == subject_genome:
                        sys.exit('Error - same query and subject genomes at this line:\n' + line)

                    for tax_level in tax_levels:
                        query_taxon = genome_taxonomy.loc[query_genome, tax_level]
                        subject_taxon = genome_taxonomy.loc[subject_genome, tax_level]

                        # Print both regions to bedfile if taxonomy at this level differs between genomes.
                        if query_taxon != subject_taxon:

                            match_summary = '||'.join([query_contig, str(query_start), str(query_end), subject_contig, str(subject_start), str(subject_end)])

                            # Print query region.
                            out_query = [query_contig, str(query_start - 1), str(query_end), 'query_' + match_summary, str(percent_identity), '+']
                            print('\t'.join(out_query), file = out_filehandles[tax_level])

                            # Determine strand of subject match and arrange accordingly.
                            if subject_start < subject_end:
                                out_subject = [subject_contig, str(subject_start - 1), str(subject_end), 'subject_' + match_summary, str(percent_identity), '+']
                            else:
                                out_subject = [subject_contig, str(subject_end - 1), str(subject_start), 'subject_' + match_summary, str(percent_identity), '-']

                            print('\t'.join(out_subject), file = out_filehandles[tax_level])

                            break


        for tax_level in tax_levels:
            out_filehandles[tax_level].close()

if __name__ == '__main__':
    main()
