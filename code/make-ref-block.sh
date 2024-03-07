#!/bin/bash

declare -A locus_start=( ["IL21"]=122121066 ["PTPN22"]=113313811 ["IFT172"]=26989805)
declare -A locus_end=( ["IL21"]=123121066 ["PTPN22"]=114313811 ["IFT172"]=27989805)
declare -A locus_chr=( ["IL21"]=4 ["PTPN22"]=1 ["IFT172"]=2)

data_dir="/home/jp2045/rds/rds-wallace-share-rU5KyGZrBOo/Data/reference"
sample_file="$data_dir/1000GP_Phase3/sparse_basis/EUR.sample"

out_file="data/$(basename ${snakemake_output[0]} .impute.legend)"
locus=$(basename $out_file .vcf.gz)

vcf_file="$data_dir/1000GP_Phase3_GRCh38/ALL.chr${locus_chr[$locus]}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"

zcat ${vcf_file} | \
  sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | \
  vcftools --gzvcf - --IMPUTE --out ${out_file} \
  --chr "chr${locus_chr[$locus]}" \
  --from-bp ${locus_start[$locus]} --to-bp ${locus_end[$locus]} \
  --remove-indels --remove-filtered-all --keep ${sample_file} \
  --maf 0.01 --max-alleles 2
