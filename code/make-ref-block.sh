#!/bin/bash

declare -A gene_tss=( ["IL21"]=122621066 ["PTPN22"]=113871759 ["IFT172"]=27489805)
declare -A gene_chr=( ["IL21"]=4 ["PTPN22"]=1 ["IFT172"]=2)

data_dir="/home/jp2045/rds/rds-wallace-share-rU5KyGZrBOo/Data/reference"
sample_file="$data_dir/1000GP_Phase3/sparse_basis/EUR.sample"

out_file="data/$(basename ${snakemake_output[0]} .impute.legend)"
gene=$(basename $out_file .vcf.gz)

width=500000
region_start=$((${gene_tss[$gene]} - $width))
region_end=$((${gene_tss[$gene]} + $width))

vcf_file="$data_dir/1000GP_Phase3_GRCh38/ALL.chr${gene_chr[$gene]}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"

echo $region_start

echo $region_end

zcat ${vcf_file} | \
  sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | \
  vcftools --gzvcf - --IMPUTE --out ${out_file} \
  --chr "chr${gene_chr[$gene]}" \
  --from-bp ${region_start[$gene]} --to-bp ${region_end[$gene]} \
  --remove-indels --remove-filtered-all --keep ${sample_file} \
  --maf 0.01 --max-alleles 2
