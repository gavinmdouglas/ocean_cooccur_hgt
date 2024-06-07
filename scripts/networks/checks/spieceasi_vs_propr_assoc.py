import gzip
import sys

# Get simplified table of association values for SPIEC-EASI and propr networks, for a quick comparison.

propr_assoc = dict()
spieceasi_assoc = dict()

with gzip.open("/mfs/gdouglas/projects/water_mags/coverm/network_working/metaG_propr_rpkm.tsv.gz", "rt") as propr_fh:
    propr_fh.readline()
    for propr_line in propr_fh:
        propr_line = propr_line.strip().split("\t")
        propr_id = "_".join(sorted(propr_line[:2]))
        propr_assoc[propr_id] = propr_line[2]

with gzip.open("/mfs/gdouglas/projects/water_mags/coverm/network_working/metaG_spieceasi_rpkm.tsv.gz", "rt") as spieceasi_fh:
    spieceasi_fh.readline()
    for spieceasi_line in spieceasi_fh:
        spieceasi_line = spieceasi_line.strip().split("\t")
        spieceasi_id = "_".join(sorted(spieceasi_line[:2]))
        spieceasi_assoc[spieceasi_id] = spieceasi_line[2]

print("taxa_combo\tpropr_assoc\tspieceasi_assoc")
for assoc_id in propr_assoc:
    if assoc_id in spieceasi_assoc:
        print(assoc_id, propr_assoc[assoc_id], spieceasi_assoc[assoc_id])

propr_assoc_id = set(list(propr_assoc.keys()))
spieceasi_assoc_id = set(list(spieceasi_assoc.keys()))

print("Propr specific:", file=sys.stderr)
print(len(propr_assoc_id - spieceasi_assoc_id), file=sys.stderr)

print("SPIEC-EASI specific:", file=sys.stderr)
print(len(spieceasi_assoc_id - propr_assoc_id), file=sys.stderr)

# Intersected:
print("Intersected:", file=sys.stderr)
print(len(propr_assoc_id & spieceasi_assoc_id), file=sys.stderr)
