#!/bin/bash
# compare the reads coverage and report intervals where resistant reads are mapped
# but susceptible reads are not mapped
res="WRC17"
samples="Hitome RIL33 RIL36 RIL48 RIL52 RIL54 RIL71"
# Change this to a bigger number to merge intervals farther away from each other
flanking=1000
# Change this number to report only intervals >= threshold
threshold=1000
for sample in $samples; do
	bedtools merge -i ${sample}.RGA.bam -d ${flanking} >> susceptible.bed
done
sort -k1,1 -k2,2n susceptible.bed > foo
mv foo susceptible.bed
bedtools merge -i susceptible.bed -d ${flanking} > susceptible_merged.bed
bedtools merge -i ${res}.RGA.bam -d ${flanking} > resistant.bed
#bedtools intersect -f 0.9 -r -v -a resistant.bed -b susceptible_merged.bed > candidates.bed
bedtools subtract -f 0.9 -r -a resistant.bed -b susceptible_merged.bed | awk '{if ($3 - $2 >= '$threshold'){ print;}}' > candidates.bed
