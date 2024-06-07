#!/usr/bin/python3

import argparse
import sys
import pandas as pd
import numpy as np
import gzip

def determine_best_hit(hits):

    highest_bitscore = -1

    for hit in hits:

        if hit['bitscore'] > highest_bitscore:

            highest_bitscore = hit['bitscore']

            list_str = []
            for index, value in hit.items():
                list_str.append(str(value))

            best_hit = list_str

    return(best_hit)


def process_blast_hit(blast_table, scaffold_column, start_column, end_column):
    '''
    Workflow for parsing BLAST hits and identifying clusters of overlapping hits.
    Meant to be used to parse either query or subject coordinates.
    '''

    num_unique_clusters = 0
    multi_cluster_sizes = []

    passed_hits = []

    for unique_contig in sorted(blast_table[scaffold_column].unique()):

        blast_table_subset = blast_table.loc[blast_table[scaffold_column] == unique_contig, :]

        # Then sort table based on subject name and hit lower and upper coordinates (in ascending order).
        blast_table_subset = blast_table_subset.sort_values([start_column, end_column],
                                                                ascending = True)
        previous_higher = -1
        cluster_hits = []

        for i, row in blast_table_subset.iterrows():
     
            current_scaffold = row[scaffold_column]
            current_lower = row[start_column]

            if current_scaffold != unique_contig:
                sys.exit('Error - Scaffold ID does not match current contig!')

            if current_lower > previous_higher:

                cluster_size = len(cluster_hits)

                if cluster_size > 0:
                    cluster_best_hit = determine_best_hit(cluster_hits)
                    passed_hits.append(cluster_best_hit)
                    num_unique_clusters += 1

                    if cluster_size > 1:
                        multi_cluster_sizes.append(cluster_size)

                cluster_hits = []

            cluster_hits.append(row)

            if row[end_column] > previous_higher:
                previous_higher = row[end_column]

        # Print final cluster.
        cluster_size = len(cluster_hits)
        if cluster_size > 0:
            cluster_best_hit = determine_best_hit(cluster_hits)
            passed_hits.append(cluster_best_hit)
            num_unique_clusters += 1

            if (cluster_size > 1):
                multi_cluster_sizes.append(cluster_size)

    return(passed_hits, num_unique_clusters, multi_cluster_sizes)


def main():

    parser = argparse.ArgumentParser(

        description=
        'Keep only a single BLAST hit per cluster of overlapping BLAST hits. '
        'Does this by first comparing hits with overlapping subject coordinates, '
        'choosing a single best hit, and then performing the same steps based on '
        'overlapping query coordinates. '
        'Very importantly, this script assumes that the BLAST output is in this order: '
        'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore, '
        'and these columns can optionally be at the end: qlen and slen.',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--blast_results", metavar="FRAGMENT_BED", type=str,
                        help="BLAST results ", required=True)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Output BLAST file", required=True)

    parser.add_argument('--no_log_header', default=False, action='store_true',
                        help='Flag to suppress header in log output. Useful if you want to combine the log of many outputs together.')


    args = parser.parse_args()

    # Check if blast results file is empty or not. If empty, then end script early and explain why.
    # Do not exit with error -- gracefully exit.
    if sum(1 for line in gzip.open(args.blast_results, 'rt')) == 0:
        print('BLAST results file is empty. Exiting gracefully.', file = sys.stderr)
    else:
        blast_results = pd.read_table(filepath_or_buffer = args.blast_results,
                                    sep = "\t",
                                    header = None)

        if len(blast_results.columns) == 14:
            blast_results.columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".split(" ")
        elif len(blast_results.columns) == 12:
            blast_results.columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")
        else:
            sys.exit("Stopping - number of columns is different from expected.")

        # Add in additional column of lower subject position, as that's the important point when sorting.
        blast_results['subject_lower'] = -99999
        blast_results['subject_higher'] = -99999
        for i, row in blast_results.iterrows():
            if row['sstart'] < row['send']:
                blast_results.loc[i, 'subject_lower'] = row['sstart']
                blast_results.loc[i, 'subject_higher'] = row['send']
            elif row['sstart'] > row['send']: 
                blast_results.loc[i, 'subject_lower'] = row['send']
                blast_results.loc[i, 'subject_higher'] = row['sstart']
            else:
                sys.exit('Error comparing subject coordinates.')

        # Loop through, identify clusters of overlapping hits, and retain hits to then be checked usng the same approach for query coordinates.
        subject_processed_out = process_blast_hit(blast_table = blast_results,
                                                scaffold_column = 'sseqid',
                                                start_column = 'subject_lower',
                                                end_column = 'subject_higher')
        
        subject_passed_hits = subject_processed_out[0]
        subject_num_unique_clusters = subject_processed_out[1]
        subject_multi_cluster_sizes = subject_processed_out[2]

        # Then convert list of lists to pandas dataframe, with same original header.
        # And with correct data types per column.
        subject_passed_hits = pd.DataFrame(subject_passed_hits)
        subject_passed_hits.columns = blast_results.columns
        for column_id in subject_passed_hits.columns:
            if column_id not in ['qseqid', 'sseqid']:
                subject_passed_hits[column_id] = pd.to_numeric(subject_passed_hits[column_id])
            else:
                subject_passed_hits[column_id] = subject_passed_hits[column_id].astype(str)

        # Check if there are any cases of the query start coordinate being higher than the query end coordinate.
        if len(subject_passed_hits.loc[subject_passed_hits['qstart'] > subject_passed_hits['qend'], :]) > 0:
            sys.exit('Error - query start coordinate is higher than query end coordinate.')

        # Then repeat the same process, but for query coordinates.
        query_processed_out = process_blast_hit(blast_table = subject_passed_hits,
                                                scaffold_column = 'qseqid',
                                                start_column = 'qstart',
                                                end_column = 'qend')
        
        query_passed_hits = query_processed_out[0]
        query_num_unique_clusters = query_processed_out[1]
        query_multi_cluster_sizes = query_processed_out[2]

        with open(args.output, 'w') as out_fh:
            for passed_hit in query_passed_hits:
                print('\t'.join([str(x) for x in passed_hit]), file = out_fh)

        # Print out log information.
        log_colnames = ['infile',
                        'number_hits',
                        'number_unique_subject_clusters',
                        'number_multi_subject_clusters',
                        'mean_multi_subject_cluster_size',
                        'sd_multi_subject_cluster_size',
                        'number_unique_query_clusters_post_subject_dedup',
                        'number_multi_query_clusters_post_subject_dedup',
                        'mean_multi_query_cluster_size_post_subject_dedup',
                        'sd_multi_query_cluster_size_post_subject_dedup']

        if len(subject_multi_cluster_sizes) > 0:
            mean_subject_multi_cluster_size = np.mean(subject_multi_cluster_sizes)
            sd_subject_multi_cluster_size = np.std(subject_multi_cluster_sizes)
        else:
            mean_subject_multi_cluster_size = float("nan")
            sd_subject_multi_cluster_size = float("nan")

        # Same computations, but for query coordinates.
        if len(query_multi_cluster_sizes) > 0:
            mean_query_multi_cluster_size = np.mean(query_multi_cluster_sizes)
            sd_query_multi_cluster_size = np.std(query_multi_cluster_sizes)
        else:
            mean_query_multi_cluster_size = float("nan")
            sd_query_multi_cluster_size = float("nan")

        if not args.no_log_header:
            print('\t'.join(log_colnames), file = sys.stderr)

        print('\t'.join([args.blast_results,
                        str(len(blast_results.index)),
            
                        str(subject_num_unique_clusters),
                        str(len(subject_multi_cluster_sizes)),
                        str(mean_subject_multi_cluster_size),
                        str(sd_subject_multi_cluster_size),

                        str(query_num_unique_clusters),
                        str(len(query_multi_cluster_sizes)),
                        str(mean_query_multi_cluster_size),
                        str(sd_query_multi_cluster_size)]),
            file = sys.stderr)


if __name__ == '__main__':
    main()
