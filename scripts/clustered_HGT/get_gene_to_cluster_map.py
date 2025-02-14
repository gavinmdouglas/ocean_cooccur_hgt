#!/usr/bin/python3

import sys

def main():

    # Parse CD-HIT output file of cluster definitions, and get simple table of gene to cluster ID.''',

    output_path = '/mfs/gdouglas/projects/ocean_mags/clusters/gene_to_cluster.tsv'
    cluster_path = '/mfs/gdouglas/projects/ocean_mags/clusters/gene-catalog.output.clstr'

    current_cluster = None
    past_ids = set()
    with open(output_path, 'w') as out_fh:
        out_fh.write('gene_id\tcluster_id\n')
        with open(cluster_path, 'r') as cluster_fh:
            for line in cluster_fh:

                if line[0] == '>':
                    current_cluster = 'c' + line.split()[1]
                else:
                    seq_id = line.split()[2][1:].split('...')[0]
                    if seq_id in past_ids:
                        sys.exit('Sequence ' + seq_id + ' appears in multiple clusters.')
                    else:
                        past_ids.add(seq_id)
                    out_fh.write(seq_id + '\t' + current_cluster + '\n')


if __name__ == '__main__':
    main()
