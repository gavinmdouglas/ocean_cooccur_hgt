from collections import defaultdict
import gzip

gene_to_cluster = {}
cluster_to_genes = dict()
gene_to_cluster = dict()
with open('/mfs/gdouglas/projects/ocean_mags/clusters/gene-catalog.output.clstr', 'r') as cluster_fh:
    for line in cluster_fh:
        if line.startswith('>Cluster'):
            cluster_id = line.strip().split()[-1]
            cluster_to_genes[cluster_id] = set()
        else:
            gene_id = line.split()[2][1:-3]
            cluster_to_genes[cluster_id].add(gene_id)
            gene_to_cluster[gene_id] = cluster_id

gene_to_rep = {}
rep_to_genes = defaultdict(set)
with gzip.open('/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog-membership.tsv.gz', 'rt') as rep_map_fh:
    header = next(rep_map_fh)
    for line in rep_map_fh:
        gene_id, rep_id = line.strip().split()
        gene_to_rep[gene_id] = rep_id
        rep_to_genes[rep_id].add(gene_id)


for rep in rep_to_genes.keys():
    rep_genes = rep_to_genes[rep]
    rep_clusters = set()
    missing = 0
    for gene in rep_genes:
        if gene in gene_to_cluster:
            rep_clusters.add(gene_to_cluster[gene])
        else:
            missing += 1
    print(rep, len(rep_genes), len(rep_clusters), missing)
