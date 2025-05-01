#!/bin/bash -i 
#SBATCH --time 2:00:00 
#SBATCH --output=logs/align_parts.%j.%a.log
#SBATCH --mem 16G

HOSTNAME=`hostname`
echo "starting (${HOSTNAME}) ${SPECIES}"

TEMPDIR=`mktemp -d`

python extract.py ../../results/panther_genes_for_alignment.tsv ${FAM} ./db > ${TEMPDIR}/fam.fa
mafft --auto ${TEMPDIR}/fam.fa > ${TEMPDIR}/msa.fa
python compute_pident.py ${TEMPDIR}/msa.fa > ~/scratch/align-res/${FAM}.tsv

rm -r ${TEMPDIR}

echo "done"
