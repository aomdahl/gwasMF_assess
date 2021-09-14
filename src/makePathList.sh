#!/bin/bash
set -e
#LIST="../../gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv"
LIST=$1
while read p; do
ls /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/${p}*.both_sexes.tsv
done <$LIST
