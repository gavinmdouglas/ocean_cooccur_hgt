#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict

# Read in FASTA file.
fasta = defaultdict(str)

with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/download/representatives.aquatic.contigs.fasta.gz', 'rt') as fasta_fh:
    for line in fasta_fh:
        if line.startswith('>'):
            sequence_id = line[1:].strip()
            contig_subset = sequence_id.split(' ')[0]
        else:
            fasta[contig_subset] += line.strip()

print("Done reading in FASTA.", file=sys.stderr)

lengths = []
print("contig\tlength")
for contig in fasta:
    print(contig + '\t' + str(len(fasta[contig])))
    lengths.append(len(fasta[contig]))

# Min length:
print("Min length: " + str(min(lengths)), file=sys.stderr)