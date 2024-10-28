#!/usr/bin/env bash
# Run the filtering on the id mapping from UniProt, to identify mappings between bgee-panther.
# - bgee mostly uses ensembl gene ids, but also wormbase, etc.

URL=https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

TMPDIR=`mktemp -d`
wget -O ${TMPDIR}/idmapping.dat.gz ${URL}

PANTHER_UNIPROT_FN=../results/pthr18_all_uniprot_ids.txt
BGEE_GENEID_FN=../results/bgee_15.1_all_gene_ids.txt

OUT_FN=../results/pthr-bgee-map.tsv.gz
# run the filtering
python filter_idmapping.py ${TMPDIR}/idmapping.dat.gz ${PANTHER_UNIPROT_FN} ${BGEE_GENEID_FN} | gzip > ${OUT_FN}

# cleanup
rm -r ${TMPDIR}
