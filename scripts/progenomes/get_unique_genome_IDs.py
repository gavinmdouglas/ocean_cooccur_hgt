import gzip

genomes = set()

with gzip.open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/download/representatives.aquatic.contigs.fasta.gz', 'rt') as fasta_fh:
    for line in fasta_fh:
        if line.startswith('>'):
            sequence_id = line[1:].strip()
            contig_subset = sequence_id.split(' ')[0]
            contig_subset_split = contig_subset.split('.')
            genome_id = 'g' + contig_subset_split[0] + '_' + contig_subset_split[1]
            genomes.add(genome_id)

for genome in sorted(list(genomes)):
    print(genome)
