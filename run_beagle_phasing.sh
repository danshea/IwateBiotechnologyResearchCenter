#!/bin/bash

BEAGLE_DIR=/dsu0/LocalSoftware/beagle_v5.0
DATA_DIR=SNPs_between_Hitomebore_and_founder
SUFFIX=_imputed_SNP_filtered.vcf.gz

[[ ! -d beagle_output ]] && mkdir -p beagle_output

for datafile in ${DATA_DIR}/*${SUFFIX}; do
	sample_name=$(basename ${datafile} ${SUFFIX})
	founder=$(grep ${sample_name} ${DATA_DIR}/samples.txt | cut -f2)
	echo "$(date +%Y%m%D-%H%M%S): Analyzing ${founder}" &>> beagle_run_$(date +%Y%m%D-%H%M%S).out
	[[ ! -d beagle_output/${sample_name}_${founder} ]] && mkdir -p beagle_output/${sample_name}_${founder}
	${BEAGLE_DIR}/run_beagle.sh gt=${datafile} out=beagle_output/${sample_name}_${founder}/${sample_name} ne=3000
done
