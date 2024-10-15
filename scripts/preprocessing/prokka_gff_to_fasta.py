#!/usr/bin/python3

import argparse
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions import write_fasta, reverse_complement

def main():

    parser = argparse.ArgumentParser(

        description='''
Parse Prokka GFF file (with scaffold FASTAs at bottom), and output a FASTA of all CDS sequences for this genome.
''',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-g', '--gff',
                        metavar="GFF", type=str,
                        help="Path to GFF file.",
                        required=True)

    parser.add_argument('-o', '--outfile',
                        metavar="OUTFILE", type=str,
                        help="Path to output file.",
                        required=True)

    args = parser.parse_args()

    # Read through once to read in all scaffolds at bottom.
    fasta_line = False
    scaffolds = {}
    exp_scaffolds_lengths = dict()
    with open(args.gff, 'r') as gff_fh:
        for line in gff_fh:
            if line.startswith('##sequence-region'):
                line_split = line.split()
                exp_scaffolds_lengths[line_split[1]] = int(line_split[3])
            elif line.startswith('##FASTA'):
                fasta_line = True
                continue
            elif fasta_line:
                if line.startswith('>'):
                    scaffold_id = line.split()[0][1:]
                    scaffolds[scaffold_id] = ''
                else:
                    scaffolds[scaffold_id] += line.rstrip()

    # Sanity check that all expected scaffolds present, and that lengths match.
    for scaffold_id in sorted(exp_scaffolds_lengths.keys()):
        exp_length = exp_scaffolds_lengths[scaffold_id]
        if scaffold_id not in scaffolds.keys():
            sys.exit('Scaffold ' + scaffold_id + ' is missing from FASTA section of GFF file.')
        if len(scaffolds[scaffold_id]) != exp_length:
            sys.exit('Scaffold ' + scaffold_id + ' has length ' + str(len(scaffolds[scaffold_id])) + ' but expected length ' + str(exp_length) + '.')

    # Read through again to extract CDS sequences.
    cds = {}
    with open(args.gff, 'r') as gff_fh:
        for line in gff_fh:
            line_split = line.rstrip().split('\t')
            if len(line_split) < 7:
                continue
            elif line_split[2] != 'CDS':
                continue

            scaffold_id = line_split[0]
            start = int(line_split[3])
            end = int(line_split[4])
            strand = line_split[6]
            cds_seq = scaffolds[scaffold_id][start-1:end]

            if strand == '-':
                cds_seq = reverse_complement(cds_seq)
            elif strand != '+':
                sys.exit('Strand not + or -.')

            gene_id = line_split[8].split(';')[0].split('=')[1]

            if gene_id in cds.keys():
                sys.exit('Gene ID ' + gene_id + ' is duplicated.')
            else:
                cds[gene_id] = cds_seq

    # Write out CDS sequences.
    write_fasta(cds, args.outfile)

if __name__ == '__main__':
    main()
