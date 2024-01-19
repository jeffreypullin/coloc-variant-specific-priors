#!/bin/bash

source ~/.bashrc
module load R/4.2.2
micromamba activate snakemake

local_flag=''

print_usage() {
  printf "Usage: -l to run locally otherwise submits to server\n"
}

while getopts 'l' flag; do
  case "${flag}" in
    l) local_flag='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [[ "$local_flag" = "true" ]]
then
  snakemake -s Snakefile --profile none --workflow-profile none --cores 1
else
  snakemake -s Snakefile --profile none
fi


