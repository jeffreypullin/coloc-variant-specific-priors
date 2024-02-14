#!/bin/bash

source ~/.bashrc
micromamba activate vcftools

data_dir="/home/jp2045/rds/rds-wallace-share-rU5KyGZrBOo/Data/reference"
vcf_file="$data_dir/1000GP_Phase3_GRCh38/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
sample_file="$data_dir/1000GP_Phase3/sparse_basis/EUR.sample"

out_file="data/$(basename ${snakemake_output[0]} .impute.legend)"

zcat ${vcf_file} | \
  sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | \
  vcftools --gzvcf - --IMPUTE --out ${out_file} \
  --chr chr1 --from-bp 113671759 --to-bp  114071759 \
  --remove-indels --remove-filtered-all --keep ${sample_file} \
  --maf 0.01 --max-alleles 2