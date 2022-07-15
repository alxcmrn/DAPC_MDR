#!/usr/bin/env bash

usage="\

Filter SNP data sets for DAPC missing data study
Generates vcfs with 5%, 10% and 20% missing data
Author: Alex Cameron

Usage:
  ./$(basename "$0") data.vcf Output_prefix
  
  (1) Vcf file to be filtered
  (2) Prefix to be appended to 3 output vcfs

"

if { [[ $# -ne 2 ]]; };then
    echo -n "$usage">&1
    exit 1
fi

inputvcf=$1
outputprefix=$2

## max-missing = 50%
## minor allele count = 3
## min quality = 30
## min depth = 5

if [[ ${inputvcf} = *.gz ]]; then
vcftools --gzvcf ${inputvcf} --max-missing 0.5 --mac 3 --minQ 30 --minDP 5 --recode --recode-INFO-all --out temp1
else
vcftools --vcf ${inputvcf} --max-missing 0.5 --mac 3 --minQ 30 --minDP 5 --recode --recode-INFO-all --out temp1
fi

## remove individuals with greater than 30% missing data
vcftools --vcf temp1.recode.vcf --missing-indv
mawk '$5 >= 0.30' out.imiss | cut -f1 | sed '1,1d'> lowDP.indv
vcftools --vcf temp1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out temp2

## use vcflib to decompose multiallelic states to single variants; remove indels and isolate diallelic SNPs
vcfallelicprimitives -k -g temp2.recode.vcf | sed 's:\.|\.:\.\/\.:g' > temp2.5.vcf
vcftools --vcf temp2.5.vcf --max-alleles 2 --min-alleles 2 --remove-indels --recode --recode-INFO-all --out temp3

## minor allele freq of 5%
vcftools --vcf temp3.recode.vcf --maf 0.05 --recode --recode-INFO-all --out temp4

## remove loci with 2X mean depth 
bgzip temp4.recode.vcf && tabix temp4.recode.vcf.gz
bcftools query temp4.recode.vcf.gz -f  '%DP\n' > depth
depth_cutoff=`awk -F'\t' '{sum+=$1}END{print 2*(sum/NR)}' depth`
bcftools filter -e "INFO/DP > ${depth_cutoff}" -o temp5.vcf temp4.recode.vcf.gz

## thin data keeping one SNP per 1000 bp
vcftools --vcf temp5.vcf --thin 1000 --recode --recode-INFO-all --out temp6

## generate data sets with 5% 20% and 30% missing data
echo "Generating data sets with variable amounts of missing data..."
vcftools --vcf temp6.recode.vcf --max-missing 0.95 --recode --recode-INFO-all --out ${outputprefix}_05per
vcftools --vcf temp6.recode.vcf --max-missing 0.90 --recode --recode-INFO-all --out ${outputprefix}_10per
vcftools --vcf temp6.recode.vcf --max-missing 0.80 --recode --recode-INFO-all --out ${outputprefix}_20per

rm temp*.log temp*.vcf temp4*

## gather SNP counts

if [[ -f ${outputprefix}_SNPcounts.txt ]]; then
	rm ${outputprefix}_SNPcounts.txt
fi

for file in $(ls *per.recode.vcf); do mawk '!/#/' $file | wc -l | paste - <(echo "$file") >> ${outputprefix}_SNPcounts.txt; done

bgzip ${outputprefix}_05per.recode.vcf
bgzip ${outputprefix}_20per.recode.vcf
bgzip ${outputprefix}_10per.recode.vcf

echo "Done"
