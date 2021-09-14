#!/bin/bash
set -e
REF="/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/ldsc_reference/"
FILES=$1
OUTDIR=$2 #relative path is fine
CURR_DIR=`pwd`
NFILE=$3 #specify the number of commands per output file.
COUNTER=0
FNAME=0
while read p; do
COND=$[ $COUNTER % $NFILE ]
if [ $COND -eq "0" ]
then
    FNAME=$((FNAME+1))
    echo "ml python/3.7-anaconda" > ${FNAME}.run.sh
    echo "cd $REF" >> ${FNAME}.run.sh
fi
on=`basename ${p} | cut -f 1 -d "."`
echo "python2 /work-zfs/abattle4/ashton/genomics_course_2020/project_2/ldsc/ldsc.py \
        --h2-cts ${p} \
        --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
        --out ${CURR_DIR}/${OUTDIR}/${on} \
        --ref-ld-chr-cts Multi_tissue_chromatin.ldcts \
        --w-ld-chr weights_hm3_no_hla/weights." >> ${FNAME}.run.sh
COUNTER=$((COUNTER+1))
done < $FILES
