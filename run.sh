#!/bin/bash

source ~/.bashrc
module load R/4.2.0
micromamba activate snakemake

snakemake -s Snakefile.py --profile none
