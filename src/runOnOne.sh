set -e 
ml gcc/5.5.0
ml R
cd /work-zfs/abattle4/ashton/snp_networks/scratch/ldsc_all_traits
    T="/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/"
fact=$1    
mkdir -p results/${fact}
  Rscript src/factAssessment.R --factors ${T}/factorization_data/${fact}.factors.txt  \
    --output ./results/${fact} \
    --ldsc_reference ldsc_results/seed2_thres0.9_h2-0.1/ \
    --ldsc_dir ${T}/results/${fact}/ldsc_enrichment_Multi_tissue_chromatin/ \
    --trait.ids /work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv \
    --trait.ids /work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv
