configfile: "config.yaml"

chrs = [x for x in range(1, 23)]

rule all: 
  input: 
    # Main text figures.
    "output/figures/prior-plot.pdf",
    "output/figures/simulation-plot.pdf",
    "output/figures/gwas-eqtl-overall-impact-plot.pdf",
    "output/figures/nek6-psmb7-example-plot.pdf",
    "output/figures/pqtl-eqtl-perf-both-plot.pdf",
    "output/figures/ukbb-gwas-eqtl-plot.pdf",
    # Supplementary figures.
    "output/figures/eqtl-dist-plot.pdf",
    "output/figures/onek1k-plot.pdf",
    "output/figures/otg-probs-plot.pdf",
    "output/figures/pqtl-eqtl-prior-effect-plot.pdf",
    "output/figures/benchmark-plot.pdf",
    # Supplementary spreadsheet.
    "output/tables/gwas-eqtl-results.xlsx",
    expand("data/output/ukbb-gwas-eqtl-coloc-abf-{icd_code_eqtl_id}-{chr}.rds", icd_code_eqtl_id = config["ukbb_gwas_eqtl_coloc_ids"], chr = chrs) 

# Download metadata.

rule download_metadata:
  output: 
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    pqtl_metadata_file = "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz",
    eqtl_catalogue_metadata = "data/eqtl-catalogue/eqtl-catalogue-metadata.tsv"
  localrule: True
  shell:
    """
    wget -O {output.eqtl_metadata_file} https://zenodo.org/record/7808390/files/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz
    wget -O {output.pqtl_metadata_file} https://zenodo.org/record/7808390/files/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz
    wget -O {output.eqtl_catalogue_metadata} https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv
    """ 

rule download_tss_data: 
  output: hg19_tss_data_path = "data/tss-data/hg19-tss-data.rds",
  script: "code/download-tss-data.R"

# Download and process eQTL data.
 
rule download_eqtlgen_data:
  output: "data/eqtlgen.txt.gz",
  shell:
    """
    wget -O {output} https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
    """

rule process_eqtlgen_data:
  input: 
    eqtlgen_path =  "data/eqtlgen.txt.gz",
    hg19_tss_data_path = "data/tss-data/hg19-tss-data.rds",
  output: processed_data_path = "data/processed-data/eqtlgen.rds"
  script: "code/process-data/process-eqtlgen-data.R"
    
rule download_eqtl_catalogue_files: 
  input: metadata_file = "data/eqtl-catalogue/eqtl-catalogue-metadata.tsv"
  output: 
    dataset_file = "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz",
    tabix_file =  "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz.tbi",
    permuted_file = "data/eqtl-catalogue/sumstats/{dataset_id}.permuted.tsv.gz",
    cred_set_file = "data/eqtl-catalogue/susie/{dataset_id}.credible_sets.tsv.gz",
    lbf_file = "data/eqtl-catalogue/susie/{dataset_id}.lbf_variable.txt.gz"
  script: "code/download-eqtl-catalogue.R"
  
rule filter_eqtl_catalogue_files:
  input: "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz"
  output: "data/eqtl-catalogue/processed-sumstats/{dataset_id}.cc.tsv"
  shell:
    """
    zcat {input} | awk 'NR == 1 {{ print }} NR != 1 {{ if ($9 <= 5E-8) {{ print }} }} ' > {output}
    """

rule run_tabix_susie_lbf_eqtl_catalogue_files:
  input: "data/eqtl-catalogue/susie/{dataset_id}.lbf_variable.txt.gz",
  output: "data/eqtl-catalogue/susie/{dataset_id}.lbf_variable.txt.gz.tbi",
  resources: 
    mem_mb = 30000,
    time_min = 120 
  shell:
    """
    mkdir -p data/tmp
    cp {input} data/tmp/{wildcards.dataset_id}.gz
    zcat data/tmp/{wildcards.dataset_id}.gz | awk -F'\\t' 'NR==1 {{print $0; next}} {{print $0 | "sort -S20G -k4,4n -k5,5n"}}' | bgzip > data/tmp/{wildcards.dataset_id}-sorted.gz
    tabix -s 4 -b 5 -e 5 -S 1 -f data/tmp/{wildcards.dataset_id}-sorted.gz
    mv data/tmp/{wildcards.dataset_id}-sorted.gz {input}
    mv data/tmp/{wildcards.dataset_id}-sorted.gz.tbi {output}
    rm data/tmp/{wildcards.dataset_id}.*
    """

rule process_eqtl_catalogue_data:
  input:
    eqtl_catalogue_metadata = "data/eqtl-catalogue/eqtl-catalogue-metadata.tsv",
    tss_data_path = "data/tss-data/hg38-tss-data.rds",
    eqtl_catalogue_paths = expand("data/eqtl-catalogue/processed-sumstats/{dataset_id}.cc.tsv",
                                  dataset_id = config["eqtl_catalogue_dataset_ids"]),
  output:
    processed_data_path = "data/processed-data/eqtl-catalogue.rds"
  script: "code/process-data/process-eqtl-catalogue-data.R"
  
rule download_onek1k_data:
  output: "data/onek1k.tsv.gz"
  shell:
    """
    wget -O {output} https://onek1k.s3.ap-southeast-2.amazonaws.com/esnp/esnp_table.tsv.gz
    """
    
rule process_onek1k_data:
  input: 
    onek1k_path =  "data/onek1k.tsv.gz",
    tss_data_path = "data/tss-data/hg19-tss-data.rds"
  output: processed_data_path = "data/processed-data/onek1k.rds"
  script: "code/process-data/process-onek1k-data.R"

# Run PolyFun

rule download_annotation_data:
    output: "data/polyfun-precomputed-annotations.tar.gz"
    shell:
      """
      wget -O {output} "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.polyfun.tar.gz"
      gunzip polyfun-precomputed-annotations.tar.gz
      tar tvf polyfun-precomputed-annotations.tar.gz 
      """

rule download_ukbb_saige_data:
    output: "data/ukbb-saige/phenocode-{icd_code}.tsv.gz"
    shell: 
      """
      wget -O {output} "https://pheweb.org/UKB-SAIGE/download/{wildcards.icd_code}"
      """

rule download_liftover_file:
    output: "data/liftover/hg19ToHg38.over.chain.gz"
    shell:
       """
       wget -O {output} https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
       """

rule liftover_hg19tohg38_ukbb_saige_data:
    input: 
        summary_stats = "data/ukbb-saige/phenocode-{icd_code}.tsv.gz",  
        chain_file = "data/liftover/hg19ToHg38.over.chain.gz"
    output: "data/ukbb-saige/hg38-phenocode-{icd_code}.tsv.gz"
    # We download the file to avoid accidentally changing the SAIGE
    # UKBB summary statistics files, which causes PolyFun to be re-run.
    shell:
      """ 
      wget -O data/ukbb-saige/tmp-{wildcards.icd_code}.tsv.gz https://pheweb.org/UKB-SAIGE/download/{wildcards.icd_code}

      gzip -df data/ukbb-saige/tmp-{wildcards.icd_code}.tsv.gz

      awk 'BEGIN {{OFS="\t"}} NR>1 {{print "chr"$1, $2-1, $2, $0}}' data/ukbb-saige/tmp-{wildcards.icd_code}.tsv > \
       data/ukbb-saige/tmp-{wildcards.icd_code}.bed

      liftOver -bedPlus=3 -tab data/ukbb-saige/tmp-{wildcards.icd_code}.bed \
        {input.chain_file} \
        data/ukbb-saige/hg38-tmp-{wildcards.icd_code}.tsv \
        data/ukbb-saige/unmapped-{wildcards.icd_code}.bed

      awk 'BEGIN {{OFS="\t"}} {{print $4, $3, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16 }}' \
        data/ukbb-saige/hg38-tmp-{wildcards.icd_code}.tsv > \
        data/ukbb-saige/hg38-tmp-unnamed-{wildcards.icd_code}.tsv   
     
      sed "1s/.*/chrom\tpos\tref\talt\trsids\tnearest_genes\tconsequence\tpval\tbeta\tsebeta\taf\tac\ttstat/" \
        data/ukbb-saige/hg38-tmp-unnamed-{wildcards.icd_code}.tsv > \
        data/ukbb-saige/hg38-phenocode-{wildcards.icd_code}.tsv  

      rm data/ukbb-saige/hg38-tmp-{wildcards.icd_code}.tsv
      rm data/ukbb-saige/hg38-tmp-unnamed-{wildcards.icd_code}.tsv
      rm data/ukbb-saige/tmp-{wildcards.icd_code}.tsv
      rm data/ukbb-saige/tmp-{wildcards.icd_code}.bed

      gzip -f data/ukbb-saige/hg38-phenocode-{wildcards.icd_code}.tsv
      """

rule tabix_hg38_ukbb_data:
    input:  "data/ukbb-saige/hg38-phenocode-{icd_code}.tsv.gz"
    output: "data/ukbb-saige/hg38-phenocode-{icd_code}.tsv.gz.tbi"
    shell:
      """
      cp {input} data/ukbb-saige/tmp-tabix-{wildcards.icd_code}.tsv.gz
      zcat data/ukbb-saige/tmp-tabix-{wildcards.icd_code}.tsv.gz | \
        sort -k1,1n -k2,2n | \
        bgzip > data/ukbb-saige/tmp-tabix-{wildcards.icd_code}-sorted.tsv.gz
      tabix -S 1 -s 1 -b 2 -e 2 data/ukbb-saige/tmp-tabix-{wildcards.icd_code}-sorted.tsv.gz
      rm data/ukbb-saige/tmp-tabix-{wildcards.icd_code}.tsv.gz
      mv data/ukbb-saige/tmp-tabix-{wildcards.icd_code}-sorted.tsv.gz {input}
      mv data/ukbb-saige/tmp-tabix-{wildcards.icd_code}-sorted.tsv.gz.tbi {output}
      """

rule rename_ukbb_saige_data:
    input: "data/ukbb-saige/phenocode-{icd_code}.tsv.gz"
    output: "data/ukbb-saige/renamed-phenocode-{icd_code}.tsv.gz"
    shell:
      """
        gzip -cd "data/ukbb-saige/phenocode-{wildcards.icd_code}.tsv.gz" | \
            sed "1s/.*/CHROM\tPOS\tA1\tA0\tRSID\tnearest_genes\tconsequence\tP\tBETA\tSE\tMAF\tac\ttstat/" | \
            gzip > "data/ukbb-saige/renamed-phenocode-{wildcards.icd_code}.tsv.gz"
      """

rule munge_ukbb_saige_sumstats: 
    input:  "data/ukbb-saige/renamed-phenocode-{icd_code}.tsv.gz"
    output: "data/ukbb-saige/munged-phenocode-{icd_code}.parquet"
    params:
      n = lambda wildcards: {"250.2": 18945, "401": 77977, "244": 14871}[wildcards.icd_code]
    resources:
      mem_mb = 20000,
      time_min = 60
    shell: 
      """
      python3 polyfun/munge_polyfun_sumstats.py \
         --sumstats data/ukbb-saige/renamed-phenocode-{wildcards.icd_code}.tsv.gz \
         --n {params.n} \
         --out "data/ukbb-saige/munged-phenocode-{wildcards.icd_code}.parquet" \
         --min-info 0 \
         --min-maf 0.01
      """

rule run_polyfun_l2_reg_s_ldsc:
    input: "data/ukbb-saige/munged-phenocode-{icd_code}.parquet"
    output: "data/polyfun-output/step-2-{icd_code}.txt"
    resources:
        mem_mb = 120000,
        time_min = 240
    shell:
      """
      source ~/.bashrc
      micromamba activate polyfun
      mkdir -p data/polyfun-output/{wildcards.icd_code}
      python polyfun/polyfun.py \
        --compute-h2-L2 \
        --allow-missing \
        --output-prefix data/polyfun-output/{wildcards.icd_code}/polyfun-out  \
        --sumstats "data/ukbb-saige/munged-phenocode-{wildcards.icd_code}.parquet" \
        --ref-ld-chr data/baselineLF2.2.UKB/baselineLF2.2.UKB. \
        --w-ld-chr data/baselineLF2.2.UKB/weights.UKB.
      micromamba deactivate
      touch data/polyfun-output/step-2-{wildcards.icd_code}.txt
      """

rule run_polyfun_ld_scores_snp_bin:
    input: "data/polyfun-output/step-2-{icd_code}.txt"
    output: "data/polyfun-output/step-3-{icd_code}-{chr}.txt"
    resources:
        mem_mb = 30000,
        time_min = 240,
        tmpdir = "/home/jp2045/coloc-variant-specific-priors/data",
        # This limits the number of jobs submitted to 10 when 
        # snakemake is called with --resources load=10000
        # load is set to a high value so as not to affect 
        # other rules. Limiting the number of jobs is necessary 
        # to prevent running out of disk space.
        load = 1000
    shell:
      """
      mkdir -p data/polyfun-output
      source ~/.bashrc
      micromamba activate polyfun
      python polyfun/polyfun.py \
        --compute-ldscores \
        --output-prefix data/polyfun-output/{wildcards.icd_code}/polyfun-out \
        --ld-ukb \
        --chr {wildcards.chr}
      micromamba deactivate
      touch "data/polyfun-output/step-3-{wildcards.icd_code}-{wildcards.chr}.txt"
      """

rule reestimate_per_snp_heritabilities:
   input: expand("data/polyfun-output/step-3-{icd_code}-{chr}.txt", icd_code = config["ukbb_saige_icd_codes"], chr = chrs)
   output: "data/polyfun-output/step-4-{icd_code}.txt"
   resources:
     mem_mb = 30000,
     time_min = 240
   shell:
      """
      source ~/.bashrc
      micromamba activate polyfun
      python polyfun/polyfun.py \
        --compute-h2-bins \
        --allow-missing \
        --output-prefix data/polyfun-output/{wildcards.icd_code}/polyfun-out \
        --sumstats "data/ukbb-saige/munged-phenocode-{wildcards.icd_code}.parquet" \
        --w-ld-chr data/baselineLF2.2.UKB/weights.UKB.
      micromamba deactivate
      touch "data/polyfun-output/step-4-{wildcards.icd_code}.txt"
      """

rule liftover_hg19tohg38_polyfun_output:
    input: 
        summary_stats = "data/polyfun-output/{icd_code}/polyfun-out.{chr}.snpvar_constrained.gz",
        chain_file = "data/liftover/hg19ToHg38.over.chain.gz"
    output: "data/polyfun-output/{icd_code}/hg38-polyfun-out.{chr}.snpvar_constrained.gz"
    # TODO Move to external file.
    shell:
      """
      cp {input.summary_stats} data/polyfun-output/cp-{wildcards.icd_code}-{wildcards.chr}.gz

      gzip -df data/polyfun-output/cp-{wildcards.icd_code}-{wildcards.chr}.gz

      awk 'BEGIN {{OFS="\t"}} NR>1 {{print "chr"$1, $3-1, $3, $0}}' \
       data/polyfun-output/cp-{wildcards.icd_code}-{wildcards.chr} > \
       data/polyfun-output/tmp-{wildcards.icd_code}-{wildcards.chr}.bed

      liftOver -bedPlus=3 -tab data/polyfun-output/tmp-{wildcards.icd_code}-{wildcards.chr}.bed \
        {input.chain_file} \
        data/polyfun-output/tmp-hg38-raw-{wildcards.icd_code}-{wildcards.chr}.tsv \
        data/polyfun-output/unmapped-{wildcards.icd_code}-{wildcards.chr}.bed

      awk 'BEGIN {{OFS="\t"}} {{print $4, $5, $3, $7, $8, $9, $10, $11, $12 }}' \
        data/polyfun-output/tmp-hg38-raw-{wildcards.icd_code}-{wildcards.chr}.tsv > \
        data/polyfun-output/tmp-hg38-unnamed-{wildcards.icd_code}-{wildcards.chr}.tsv   
     
      sed "1s/.*/CHR\tSNP\tBP\tA1\tA2\tSNPVAR\tMAF\tN\tZ/" \
        data/polyfun-output/tmp-hg38-unnamed-{wildcards.icd_code}-{wildcards.chr}.tsv > \
        data/polyfun-output/{wildcards.icd_code}/hg38-polyfun-out.{wildcards.chr}.snpvar_constrained

      rm data/polyfun-output/cp-{wildcards.icd_code}-{wildcards.chr}
      rm data/polyfun-output/tmp-{wildcards.icd_code}-{wildcards.chr}.bed
      rm data/polyfun-output/tmp-hg38-raw-{wildcards.icd_code}-{wildcards.chr}.tsv
      rm data/polyfun-output/tmp-hg38-unnamed-{wildcards.icd_code}-{wildcards.chr}.tsv   

      gzip -f data/polyfun-output/{wildcards.icd_code}/hg38-polyfun-out.{wildcards.chr}.snpvar_constrained  
      """

# Download and create prior probabilities data.

rule compute_eqtl_densitites:
  input: 
    processed_eqtlgen_data_path = "data/processed-data/eqtlgen.rds",
    processed_onek1k_data_path = "data/processed-data/onek1k.rds"
  output: 
    eqtlgen_density_path = "output/densities/eqtlgen.rds",
    onek1k_r1_density_path = "output/densities/onek1k_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_round_2.rds",
  script: "code/compute-eqtl-densities.R"

rule download_gnochhi_data: 
  output: "data/raw-gnocchi-windows.bed" 
  shell: 
    """
    wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06045-0/MediaObjects/41586_2023_6045_MOESM4_ESM.zip
    unzip 41586_2023_6045_MOESM4_ESM.zip 
    gunzip Supplementary_Data_2.bed.gz 
    mv Supplementary_Data_2.bed {output}
    rm 41586_2023_6045_MOESM4_ESM.zip 
    rm Supplementary_Data_1.tsv  Supplementary_Data_3.bed.gz Supplementary_Data_4.tsv
    rm Supplementary_Data_5.tsv Supplementary_Data_6_ESM.txt
    """

rule process_gnocchi_data:
  input: "data/raw-gnocchi-windows.bed",
  output: "data/gnocchi-windows.bed"
  script: "code/process-gnocchi-data.R"

rule download_abc_score_data:
  output: "data/raw-abc-data.txt.gz"
  shell:
    """
    wget -O {output} https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
    """

rule process_abc_score_data:
  input: "data/raw-abc-data.txt.gz",
  output: "data/abc-data.txt.gz"
  resources: mem_mb = 10000
  script: "code/process-abc-data.R"

rule download_polyfun_data:
  output: 
    chr1_7 = "data/snpvar_meta.chr1_7.parquet",
    chr8_22 = "data/snpvar_meta.chr8_22.parquet"
  shell:
   """
   wget -O {output.chr1_7} https://github.com/omerwe/polyfun/raw/master/snpvar_meta.chr1_7.parquet
   wget -O {output.chr8_22} https://github.com/omerwe/polyfun/raw/master/snpvar_meta.chr8_22.parquet
   """

# Download Finngen GWAS data.

rule download_finngen_manifest:
  output: "data/finngen/finngen-manifest.tsv"
  shell:
   """
   wget -O {output} https://storage.googleapis.com/finngen-public-data-r10/summary_stats/R10_manifest.tsv
   """

rule download_finngen_gwas_data: 
  output: 
    sum_stat = "data/finngen/{gwas_id}.gz",
    tabix = "data/finngen/{gwas_id}.gz.tbi",
    cs = "data/finngen/{gwas_id}.SUSIE.cred.bgz",
    lbf = "data/finngen/{gwas_id}.SUSIE.snp.bgz",
  shell:
   """
   wget -O {output.sum_stat} https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_{wildcards.gwas_id}.gz
   wget -O {output.tabix} https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_{wildcards.gwas_id}.gz.tbi
   wget -O {output.cs} https://storage.googleapis.com/finngen-public-data-r10/finemap/full/susie/finngen_R10_{wildcards.gwas_id}.SUSIE.cred.bgz
   wget -O {output.lbf} https://storage.googleapis.com/finngen-public-data-r10/finemap/full/susie/finngen_R10_{wildcards.gwas_id}.SUSIE.snp.bgz
   """

rule run_tabix_finngen_lbf_files:
  input: "data/finngen/{gwas_id}.SUSIE.snp.bgz",
  output: "data/finngen/{gwas_id}.SUSIE.snp.bgz.tbi",
  shell:
    """
    tabix -S 1 -s 5 -b 6 -e 6 data/finngen/{wildcards.gwas_id}.SUSIE.snp.bgz 
    """

# Run colocalsiations.

rule run_pqtl_eqtl_coloc_abf_colocalisation:
  input: 
    # Input files.
    eqtl_data_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz",
    eqtl_index_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz.tbi",
    permutation_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.permuted.tsv.gz",
    pqtl_data_file = "data/eqtl-catalogue/sumstats/QTD000584.cc.tsv.gz",
    pqtl_index_file = "data/eqtl-catalogue/sumstats/QTD000584.cc.tsv.gz.tbi",
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    pqtl_metadata_file = "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz",
    # Prior data.
    gnocchi_data_path = "data/gnocchi-windows.bed",
    polyfun_data_1_7_path = "data/snpvar_meta.chr1_7.parquet",
    polyfun_data_8_22_path = "data/snpvar_meta.chr8_22.parquet",
    abc_score_data_path = "data/abc-data.txt.gz",
    eqtlgen_density_path = "output/densities/eqtlgen.rds",
    onek1k_r1_density_path = "output/densities/onek1k_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_round_2.rds",
  output: 
    result_file = "data/output/pqtl-eqtl-coloc-abf-{eqtl_id}-{chr}.rds"
  retries: 1
  resources: 
    mem_mb = lambda wildcards, attempt: 7000 * attempt,
    time_min = lambda wildcards, attempt: 20 * attempt ** 4
  script: "code/run-pqtl-eqtl-coloc-abf.R"

rule run_pqtl_eqtl_coloc_susie_colocalisation:
  input: 
    # Input files.
    eqtl_cs_file = "data/eqtl-catalogue/susie/{eqtl_id}.credible_sets.tsv.gz",
    eqtl_lbf_file = "data/eqtl-catalogue/susie/{eqtl_id}.lbf_variable.txt.gz",
    eqtl_lbf_index_file = "data/eqtl-catalogue/susie/{eqtl_id}.lbf_variable.txt.gz.tbi",
    pqtl_cs_file = "data/eqtl-catalogue/susie/QTD000584.credible_sets.tsv.gz",
    pqtl_lbf_file = "data/eqtl-catalogue/susie/QTD000584.lbf_variable.txt.gz",
    pqtl_lbf_index_file = "data/eqtl-catalogue/susie/QTD000584.lbf_variable.txt.gz.tbi",
    permutation_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.permuted.tsv.gz",
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    pqtl_metadata_file = "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz",
    # Prior data.
    gnocchi_data_path = "data/gnocchi-windows.bed",
    polyfun_data_1_7_path = "data/snpvar_meta.chr1_7.parquet",
    polyfun_data_8_22_path = "data/snpvar_meta.chr8_22.parquet",
    abc_score_data_path = "data/abc-data.txt.gz",
    eqtlgen_density_path = "output/densities/eqtlgen.rds",
    onek1k_r1_density_path = "output/densities/onek1k_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_round_2.rds",
  output: 
    result_file = "data/output/pqtl-eqtl-coloc-susie-{eqtl_id}-{chr}.rds"
  retries: 1
  resources: 
    mem_mb = lambda wildcards, attempt: 7000 * attempt,
    time_min = lambda wildcards, attempt: 20 * attempt ** 3
  script: "code/run-pqtl-eqtl-coloc-susie.R"

rule run_gwas_eqtl_coloc_abf_colocalisation:
  input: 
    # Input files.
    eqtl_data_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz",
    permutation_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.permuted.tsv.gz",
    eqtl_index_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz.tbi",
    gwas_data_file = "data/finngen/{gwas_id}.gz",
    gwas_index_file = "data/finngen/{gwas_id}.gz.tbi",
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    manifest_file = "data/finngen/finngen-manifest.tsv",
    # Prior data.
    gnocchi_data_path = "data/gnocchi-windows.bed",
    polyfun_data_1_7_path = "data/snpvar_meta.chr1_7.parquet",
    polyfun_data_8_22_path = "data/snpvar_meta.chr8_22.parquet",
    abc_score_data_path = "data/abc-data.txt.gz",
    eqtlgen_density_path = "output/densities/eqtlgen.rds",
    onek1k_r1_density_path = "output/densities/onek1k_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_round_2.rds",
  output: 
    coloc_results_file = "data/output/gwas-eqtl-coloc-abf-{gwas_id}-{eqtl_id}-{chr}.rds",
    finemapping_results_file = "data/output/gwas-eqtl-finemapping-{gwas_id}-{eqtl_id}-{chr}.rds"
  retries: 1
  resources: 
    mem_mb = lambda wildcards, attempt: 7000 * attempt,
    time_min = lambda wildcards, attempt: 20 * attempt ** 3
  script: "code/run-gwas-eqtl-coloc-abf.R"

rule run_gwas_eqtl_coloc_susie_colocalisation:
  input: 
    # Input files.
    eqtl_cs_file = "data/eqtl-catalogue/susie/{eqtl_id}.credible_sets.tsv.gz",
    eqtl_lbf_file = "data/eqtl-catalogue/susie/{eqtl_id}.lbf_variable.txt.gz",
    eqtl_lbf_index_file = "data/eqtl-catalogue/susie/{eqtl_id}.lbf_variable.txt.gz.tbi",
    gwas_cs_file = "data/finngen/{gwas_id}.SUSIE.cred.bgz",
    gwas_lbf_file = "data/finngen/{gwas_id}.SUSIE.snp.bgz",
    gwas_lbf_index_file = "data/finngen/{gwas_id}.SUSIE.snp.bgz.tbi",
    permutation_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.permuted.tsv.gz",
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    # Prior data.
    gnocchi_data_path = "data/gnocchi-windows.bed",
    polyfun_data_1_7_path = "data/snpvar_meta.chr1_7.parquet",
    polyfun_data_8_22_path = "data/snpvar_meta.chr8_22.parquet",
    abc_score_data_path = "data/abc-data.txt.gz",
    eqtlgen_density_path = "output/densities/eqtlgen.rds",
    onek1k_r1_density_path = "output/densities/onek1k_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_round_2.rds",
  output: 
    result_file = "data/output/gwas-eqtl-coloc-susie-{gwas_id}-{eqtl_id}-{chr}.rds"
  retries: 1
  resources: 
    mem_mb = lambda wildcards, attempt: 7000 * attempt,
    time_min = lambda wildcards, attempt: 20 * attempt ** 3
  script: "code/run-gwas-eqtl-coloc-susie.R"

rule run_ukbb_gwas_eqtl_coloc_abf_colocalisation:
  input: 
    # Input files.
    eqtl_data_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz",
    permutation_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.permuted.tsv.gz",
    eqtl_index_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz.tbi",
    gwas_data_file = "data/ukbb-saige/hg38-phenocode-{icd_code}.tsv.gz",
     gwas_index_file = "data/ukbb-saige/hg38-phenocode-{icd_code}.tsv.gz.tbi",
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    # Prior data.
    polyfun_data_1_7_path = "data/snpvar_meta.chr1_7.parquet",
    polyfun_data_8_22_path = "data/snpvar_meta.chr8_22.parquet",
    polyfun_trait_specific_data_path = "data/polyfun-output/{icd_code}/hg38-polyfun-out.{chr}.snpvar_constrained.gz",
    eqtlgen_density_path = "output/densities/eqtlgen.rds",
    onek1k_r1_density_path = "output/densities/onek1k_round_1.rds",
  output: 
    coloc_results_file = "data/output/ukbb-gwas-eqtl-coloc-abf-{icd_code}-{eqtl_id}-{chr}.rds",
    finemapping_results_file = "data/output/ukbb-gwas-eqtl-finemapping-{icd_code}-{eqtl_id}-{chr}.rds"
  retries: 1
  resources: 
    mem_mb = lambda wildcards, attempt: 7000 * attempt,
    time_min = lambda wildcards, attempt: 20 * attempt ** 3
  script: "code/run-ukbb-gwas-eqtl-coloc-abf.R"

# Run simulations.

rule make_ref_block:
  output: 
    leg_file = "data/{gene}.vcf.gz.impute.legend",
    haps_file = "data/{gene}.vcf.gz.impute.hap",
  script: "code/make-ref-block.sh"

rule run_simulations:
  input:
    leg_file = "data/{gene}.vcf.gz.impute.legend",
    haps_file = "data/{gene}.vcf.gz.impute.hap",
  output: sim_file = "data/output/sim-result-{gene}.rds"
  script: "code/run-simulation.R"

# Make figures.

rule plot_eqtl_tss_dist: 
  input: 
    eqtl_catalogue_data_file = "data/processed-data/eqtl-catalogue.rds",
    eqtlgen_data_file = "data/processed-data/eqtlgen.rds",
    onek1k_data_file = "data/processed-data/onek1k.rds",
  output: 
    dist_plot_file = "output/figures/eqtl-dist-plot.pdf",
    onek1k_plot_file = "output/figures/onek1k-plot.pdf",
  localrule: True
  script: "code/plot-eqtl-tss-dist.R"

rule plot_otg_data: 
  output: otg_plot_path = "output/figures/otg-probs-plot.pdf",
  localrule: True
  script: "code/plot-otg-data.R"

rule plot_priors:
  input: 
    autoimmune_gwas_path = "data/finngen/AUTOIMMUNE.gz",
    gnocchi_data_path = "data/gnocchi-windows.bed",
    abc_score_data_path = "data/abc-data.txt.gz",
    density_data_r1_path = "output/densities/onek1k_round_1.rds",
    density_data_r2_path = "output/densities/onek1k_round_2.rds",
    eqtlgen_density_data_path = "output/densities/eqtlgen.rds",
    snp_var_data_1_7_path = "data/snpvar_meta.chr1_7.parquet"
  output: 
    all_priors_plot_path = "output/figures/prior-plot.pdf",
    polyfun_priors_plot_path = "output/figures/ukbb-polyfun-prior-plot.pdf"
  localrule: True
  script: "code/plot-priors.R"

rule plot_simulations: 
  input: expand("data/output/sim-result-{gene}.rds", gene = config["simulation_genes"])
  output: "output/figures/simulation-plot.pdf"
  localrule: True
  script: "code/plot-simulations.R"

rule plot_pqtl_eqtl_colocalisatons:
  input: 
    coloc_abf_paths = expand(
      "data/output/pqtl-eqtl-coloc-abf-{eqtl_id}-{chr}.rds", 
      eqtl_id = config["pqtl_eqtl_coloc_dataset_ids"],
      chr = chrs
    ),
    coloc_susie_paths = expand(
      "data/output/pqtl-eqtl-coloc-susie-{eqtl_id}-{chr}.rds", 
      eqtl_id = config["pqtl_eqtl_coloc_dataset_ids"],
      chr = chrs
    ), 
    protein_metadata_path = "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"
  output: 
    pqtl_eqtl_perf_both_plot_path = "output/figures/pqtl-eqtl-perf-both-plot.pdf",
    prior_effect_plot_path = "output/figures/pqtl-eqtl-prior-effect-plot.pdf"
  localrule: True
  script: "code/plot-pqtl-eqtl-colocs.R"

rule plot_gwas_eqtl_colocalisatons:
  input: 
    coloc_abf_paths = expand(
      "data/output/gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds", 
      gwas_id_eqtl_id = config["gwas_eqtl_coloc_ids"],
      chr = chrs 
    ),
    coloc_susie_paths = expand(
      "data/output/gwas-eqtl-coloc-susie-{gwas_id_eqtl_id}-{chr}.rds", 
      gwas_id_eqtl_id = config["gwas_eqtl_coloc_ids"],
      chr = chrs
    ), 
  output: 
    overall_impact_plot_path = "output/figures/gwas-eqtl-overall-impact-plot.pdf",
    coloc_results_excel_path = "output/tables/gwas-eqtl-results.xlsx"
  localrule: True
  script: "code/plot-gwas-eqtl-colocs.R"

rule plot_locus_example: 
  input:
    coloc_abf_paths = expand(
      "data/output/gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds", 
      gwas_id_eqtl_id = config["gwas_eqtl_coloc_ids"],
      chr = chrs 
    ),
    autoimmune_gwas_path = "data/finngen/AUTOIMMUNE.gz",
    gtex_thyroid_path = "data/eqtl-catalogue/sumstats/QTD000341.cc.tsv.gz"
  output: 
    nek6_psmb7_plot_path = "output/figures/nek6-psmb7-example-plot.pdf",
  localrule: True
  script: "code/plot-locus.R"

rule plot_ukbb_gwas_eqtl_colocalisations:
  input: 
    coloc_paths = expand(
      "data/output/ukbb-gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds", 
      gwas_id_eqtl_id = config["ukbb_gwas_eqtl_coloc_ids"],
      chr = chrs 
    ),
    fm_paths = expand(
      "data/output/ukbb-gwas-eqtl-finemapping-{gwas_id_eqtl_id}-{chr}.rds",
      gwas_id_eqtl_id = config["ukbb_gwas_eqtl_coloc_ids"],
      chr = chrs 
    )
  output: ukbb_gwas_eqtl_plot_path = "output/figures/ukbb-gwas-eqtl-plot.pdf"
  localrule: True 
  script: "code/plot-ukbb-gwas-eqtl-colocs.R"

rule plot_benchmark: 
  output: benchmark_plot_path = "output/figures/benchmark-plot.pdf",
  localrule: True 
  script: "code/plot-coloc-benchmark.R"

