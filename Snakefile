configfile: "config.yaml"

chrs = [x for x in range(1, 23)]

rule all: 
  input: 
    # Main figures.
    "output/figures/prior-plot.pdf",
    "output/figures/simulation-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-abf-perf-both-plot.pdf",
    "output/figures/gwas-eqtl-coloc-abf-n-colocs-plot.pdf",
    # Supplementargy figures.
    "output/figures/eqtl-dist-plot.pdf",
    "output/figures/onek1k-plot.pdf",
    "output/figures/dataset-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-abf-perf-by-dataset-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-abf-perf-median-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-abf-n-coloc-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-abf-pph4-scatter-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-susie-perf-by-dataset-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-susie-perf-median-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-susie-pph4-scatter-plot.pdf",
    "output/figures/gwas-eqtl-coloc-abf-boostrap-scatter-plot.pdf",
    "output/figures/gwas-eqtl-coloc-abf-prob-sig-scatter-plot.pdf",
    "output/tables/gwas-eqtl-coloc-abf-results.xlsx",
    "output/tables/gwas-eqtl-coloc-susie-results.xlsx",
    "output/figures/gwas-eqtl-coloc-susie-prior-effect-plot.pdf",
    "output/figures/pqtl-eqtl-coloc-abf-perf-max-curve-plot.pdf"

# Download metadata.

rule download_metadata:
  output: 
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    pqtl_metadata_file = "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"
  shell:
    """
    wget -O {output.eqtl_metadata_file} https://zenodo.org/record/7808390/files/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz
    wget -O {output.pqtl_metadata_file} https://zenodo.org/record/7808390/files/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz
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
    
rule download_gtex_tar: 
  output: temporary("data/gtex-v8.tar")
  shell:
    """
    wget -O {output} https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
    """

rule extract_gtex_files: 
  input: "data/gtex-v8.tar"
  output: "data/gtex-v8/{tissue}.v8.signif_variant_gene_pairs.txt"
  shell:
    """
    tar -Oxf {input} GTEx_Analysis_v8_eQTL/{wildcards.tissue}.v8.signif_variant_gene_pairs.txt.gz | gunzip -c > {output}
    """
    
rule process_gtex_data: 
  input: gtex_paths = expand("data/gtex-v8/{tissue}.v8.signif_variant_gene_pairs.txt", tissue = config["gtex_tissues"])
  output: processed_data_path = "data/processed-data/gtex.rds"
  script: "code/process-data/process-gtex-data.R"

rule download_etl_catalogue_metadata: 
  output: "data/eqtl-catalogue/eqtl-catalogue-metadata.tsv"
  shell:
    """
    wget -O {output} https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv
    """

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
  localrule: True
  shell:
    """
    zcat {input} | awk 'NR == 1 {{ print }} NR != 1 {{ if ($9 <= 5E-8) {{ print }} }} ' > {output}
    """

rule run_tabix_susie_lbf_eqtl_catalogue_files:
  input: "data/eqtl-catalogue/susie/{dataset_id}.lbf_variable.txt.gz",
  output: "data/eqtl-catalogue/susie/{dataset_id}.lbf_variable.txt.gz.tbi",
  resources: 
    time_min = 120 
  shell:
    """
    cp {input} data/tmp-{wildcards.dataset_id}.gz
    zcat data/tmp-{wildcards.dataset_id}.gz | awk -F'\\t' 'BEGIN{{OFS=FS}}NR==1{{print $0; next}}{{print $0 | "sort -T data -k4,4n -k5,5n"}}' | bgzip > {input}
    tabix -s 4 -b 5 -e 5 -S 1 -f {input}
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

# Download and create prior probabilities data.

rule compute_eqtl_densitites:
  input: 
    processed_eqtlgen_data_path = "data/processed-data/eqtlgen.rds",
    processed_onek1k_data_path = "data/processed-data/onek1k.rds"
  output: 
    eqtlgen_density_path = "output/densities/eqtlgen.rds",
    onek1k_r1_density_path = "output/densities/onek1k_cd4nc_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_cd4nc_round_2.rds",
    onek1k_r3_density_path = "output/densities/onek1k_cd4nc_round_3.rds"
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
    onek1k_r1_density_path = "output/densities/onek1k_cd4nc_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_cd4nc_round_2.rds",
    onek1k_r3_density_path = "output/densities/onek1k_cd4nc_round_3.rds"
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
    onek1k_r1_density_path = "output/densities/onek1k_cd4nc_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_cd4nc_round_2.rds",
    onek1k_r3_density_path = "output/densities/onek1k_cd4nc_round_3.rds"
  output: 
    result_file = "data/output/pqtl-eqtl-coloc-susie-{eqtl_id}-{chr}.rds"
  retries: 1
  resources: 
    mem_mb = lambda wildcards, attempt: 8000 * attempt,
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
    onek1k_r1_density_path = "output/densities/onek1k_cd4nc_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_cd4nc_round_2.rds",
    onek1k_r3_density_path = "output/densities/onek1k_cd4nc_round_3.rds"
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
    onek1k_r1_density_path = "output/densities/onek1k_cd4nc_round_1.rds",
    onek1k_r2_density_path = "output/densities/onek1k_cd4nc_round_2.rds",
    onek1k_r3_density_path = "output/densities/onek1k_cd4nc_round_3.rds",
  output: 
    result_file = "data/output/gwas-eqtl-coloc-susie-{gwas_id}-{eqtl_id}-{chr}.rds"
  retries: 1
  resources: 
    mem_mb = lambda wildcards, attempt: 7000 * attempt,
    time_min = lambda wildcards, attempt: 20 * attempt ** 3
  script: "code/run-gwas-eqtl-coloc-susie.R"

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
    gtex_data_file = "data/processed-data/gtex.rds",
    eqtlgen_data_file = "data/processed-data/eqtlgen.rds",
    onek1k_data_file = "data/processed-data/onek1k.rds",
  output: 
    dist_plot_file = "output/figures/eqtl-dist-plot.pdf",
    onek1k_plot_file = "output/figures/onek1k-plot.pdf",
    dataset_plot_file = "output/figures/dataset-plot.pdf"
  script: "code/plot-eqtl-tss-dist.R"

rule plot_priors:
  input: 
    autoimmune_gwas_path = "data/finngen/AUTOIMMUNE.gz",
    gnocchi_data_path = "data/gnocchi-windows.bed",
    abc_score_data_path = "data/abc-data.txt.gz",
    density_data_r1_path = "output/densities/onek1k_cd4nc_round_1.rds",
    density_data_r2_path = "output/densities/onek1k_cd4nc_round_2.rds",
    eqtlgen_density_data_path = "output/densities/eqtlgen.rds",
    snp_var_data_1_7_path = "data/snpvar_meta.chr1_7.parquet"
  output: "output/figures/prior-plot.pdf"
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
    abf_perf_by_dataset_plot_path = "output/figures/pqtl-eqtl-coloc-abf-perf-by-dataset-plot.pdf",
    abf_perf_median_plot_path = "output/figures/pqtl-eqtl-coloc-abf-perf-median-plot.pdf",
    abf_perf_max_plot_path = "output/figures/pqtl-eqtl-coloc-abf-perf-max-plot.pdf",
    abf_n_colocs_plot_path = "output/figures/pqtl-eqtl-coloc-abf-n-coloc-plot.pdf",
    abf_pph4_scatter_plot_path = "output/figures/pqtl-eqtl-coloc-abf-pph4-scatter-plot.pdf",
    susie_perf_by_dataset_plot_path = "output/figures/pqtl-eqtl-coloc-susie-perf-by-dataset-plot.pdf",
    susie_perf_median_plot_path = "output/figures/pqtl-eqtl-coloc-susie-perf-median-plot.pdf",
    susie_perf_max_plot_path = "output/figures/pqtl-eqtl-coloc-susie-perf-max-plot.pdf",
    susie_pph4_scatter_plot_path = "output/figures/pqtl-eqtl-coloc-susie-pph4-scatter-plot.pdf",
    abf_perf_max_curve_plot_path = "output/figures/pqtl-eqtl-coloc-abf-perf-max-curve-plot.pdf",
    pqtl_eqtl_perf_plot_path = "output/figures/pqtl-eqtl-coloc-abf-perf-both-plot.pdf"
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
    abf_n_colocs_plot_path = "output/figures/gwas-eqtl-coloc-abf-n-colocs-plot.pdf",
    abf_bootstrap_scatter_plot_path = "output/figures/gwas-eqtl-coloc-abf-boostrap-scatter-plot.pdf",
    abf_prob_sig_scatter_plot_path = "output/figures/gwas-eqtl-coloc-abf-prob-sig-scatter-plot.pdf",
    abf_coloc_results_table_path = "output/tables/gwas-eqtl-coloc-abf-results.xlsx",
    susie_coloc_results_table_path = "output/tables/gwas-eqtl-coloc-susie-results.xlsx",
    susie_prior_effect_plot_path = "output/figures/gwas-eqtl-coloc-susie-prior-effect-plot.pdf",
    abf_prior_effect_plot_path = "output/figures/gwas-eqtl-coloc-abf-prior-effect-plot.pdf",
  localrule: True
  script: "code/plot-gwas-eqtl-colocs.R"
