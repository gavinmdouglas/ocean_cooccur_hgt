import sys
import gzip

old_vals = {}
with open('/mfs/gdouglas/tmp/propr_freeliv_earlier_ms_code.tsv', 'r') as OLD_FH:
    header = OLD_FH.readline()
    for line in OLD_FH:
        line = line.strip().split('\t')
        taxa_combo = ';'.join(sorted([line[-1], line[-2]]))
        if taxa_combo in old_vals:
            sys.exit('Duplicate taxa combo in old file: {}'.format(taxa_combo))
        old_vals[taxa_combo] = float(line[-4])

with gzip.open('/mfs/gdouglas/projects/ocean_mags/networks/filtersplit/metaG_propr_rpkm.freeliv.tsv.gz', 'rt') as NEW_FH:
    header = NEW_FH.readline()
    for line in NEW_FH:
        line = line.strip().split('\t')
        taxa_combo = ';'.join(sorted([line[0], line[1]]))
        if taxa_combo in old_vals:
            print(taxa_combo, line[2], old_vals[taxa_combo])

