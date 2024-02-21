#!/bin/bash

declare -A locus_start=( ["IL21"]=122529959 ["PTPN22"]=113771759 ["IRF5"]=128927032)
declare -A locus_end=( ["IL21"]=122729959 ["PTPN22"]=113971759 ["IRF5"]=129037032)
declare -A locus_chr=( ["IL21"]=4 ["PTPN22"]=1 ["IRF5"]=7)

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
