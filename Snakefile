configfile: "config.yaml"

chromosomes = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]

rule all: 
  input: 
    "data/T1D_Chiou_34012112_1-hg38.tsv.gz",
    "data/T1D_Chiou_34012112_1-hg38.tsv.gz.tbi",
    "data/processed-data/eqtl-catalogue.rds",
    "data/processed-data/gtex.rds",
    "data/processed-data/adipos-express-marginal.rds",
    "data/processed-data/eqtlgen.rds",
    "data/processed-data/onek1k.rds",
    "data/eqtlgen-sig.txt.gz",
    #"data/eqtlgen-all.txt.gz.tbi",
    "data/abc-data.txt.gz",
    "data/fauman-hyde/eqtlgen.txt",
    "data/tss-data/hg19-tss-data.rds",
    "data/gnocchi-windows.bed",
    "data/tambets-etal-supp.xlsx",
    "data/snpvar_meta.chr1_7.parquet",
    "data/snpvar_meta.chr8_22.parquet",
    "data/abc-data.txt.gz",
    expand("data/eqtl-catalogue/processed-sumstats/{dataset_id}.cc.tsv",
            dataset_id = config["eqtl_catalogue_dataset_ids"]),
    expand("output/data/sim-result-{gene}.rds", gene = config["simulation_genes"]),
    expand("data/adipos-express/processed-ab1-eur/{chromosome}.txt", chromosome = chromosomes),
    expand(
      "output/data/gwas-eqtl-coloc-{gwas_id_eqtl_id}-{chr}.rds", 
      gwas_id_eqtl_id = config["gwas_eqtl_coloc_ids"],
      chr = [x for x in range(1, 23)]
    ),
    expand("output/data/pqtl-eqtl-coloc-{eqtl_id}-{chr}.rds", 
           chr = [x for x in range(1, 23)],
           eqtl_id = config["pqtl_eqtl_coloc_dataset_ids"])

rule fetch_t1d_gwas_data:
  output: 
    data_file = "data/T1D_Chiou_34012112_1-hg38.tsv.gz",
    tabix_file = "data/T1D_Chiou_34012112_1-hg38.tsv.gz.tbi"
  shell:
    """ 
    cp /home/jp2045/rds/rds-basis-YRWZsjDGyaU/02-Processed/T1D_Chiou_34012112_1-hg38.tsv.gz /home/jp2045/coloc-estimated-eqtl-priors/data/
    mv data/T1D_Chiou_34012112_1-hg38.tsv.gz data/tmp-T1D_Chiou_34012112_1-hg38.tsv.gz  
    gunzip --force data/tmp-T1D_Chiou_34012112_1-hg38.tsv.gz
    bgzip data/tmp-T1D_Chiou_34012112_1-hg38.tsv
    tabix -s 3 -b 4 -e 4 -c S data/tmp-T1D_Chiou_34012112_1-hg38.tsv.gz
    mv data/tmp-T1D_Chiou_34012112_1-hg38.tsv.gz.tbi {output.tabix_file}
    mv data/tmp-T1D_Chiou_34012112_1-hg38.tsv.gz {output.data_file}
    """

rule download_tss_data: 
  output: hg19_tss_data_path = "data/tss-data/hg19-tss-data.rds",
          hg38_tss_data_path = "data/tss-data/hg38-tss-data.rds"
  script: "code/download-tss-data.R"
    
rule download_fauman_hyde_data: 
  output: eqtlgen = "data/fauman-hyde/eqtlgen.txt",
          ferkingstad = "data/fauman-hyde/ferkingstad.txt" 
  shell:
    """
    wget -O {output.eqtlgen} https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-022-04706-x/MediaObjects/12859_2022_4706_MOESM3_ESM.txt
    wget -O {output.ferkingstad} https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-022-04706-x/MediaObjects/12859_2022_4706_MOESM2_ESM.txt
    """

rule download_eqtlgen_data:
  output: 
    eqtlgen_sig = "data/eqtlgen-sig.txt.gz",
    eqtlgen_all = "data/eqtlgen-all.txt.gz"
  shell:
    """
    wget -O {output.eqtlgen_sig} https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
    wget -O {output.eqtlgen_all} https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
    """

rule tabix_eqtlgen:
  input: "data/eqtlgen-all.txt.gz"
  output: "data/eqtlgen-all.txt.gz.tbi"
  conda: "envs/tabix.yaml"
  shell:
    """
    cp data/eqtlgen-all.txt.gz data/tmp-eqtlgen-all.txt.gz
    gunzip data/tmp-eqtlgen-all.txt.gz
    bgzip data/tmp-eqtlgen-all.txt
    tabix -s 4 -b 3 -e 3 -c P data/tmp-eqtlgen-all.txt.gz
    mv data/tmp-eqtlgen-all.txt.gx.tbi {output}
    """
    
rule process_eqtlgen_data:
  input: eqtlgen_path =  "data/eqtlgen-sig.txt.gz"
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
    permuted_file = "data/eqtl-catalogue/sumstats/{dataset_id}.permuted.tsv.gz"
  localrule: True
  script: "code/download-eqtl-catalogue.R"
  
rule filter_eqtl_catalogue_files:
  input: "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz"
  output: "data/eqtl-catalogue/processed-sumstats/{dataset_id}.cc.tsv"
  localrule: True
  shell:
    """
    zcat {input} | awk 'NR == 1 {{ print }} NR != 1 {{ if ($9 <= 5E-8) {{ print }} }} ' > {output}
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
  
rule download_adipos_express_tars: 
  output: temporary("data/adipos-express-marginal-eur-by-chr.tar"),
  shell:
    """
    wget ftp://mohlkeanon:anon@rc-ns-ftp.its.unc.edu/marginal_byChr_EURonly.tar.gz -O - | gunzip -d > {output}
    """

rule extract_adipos_express_marginal_files:
  input: "data/adipos-express-marginal-eur-by-chr.tar",
  output: "data/adipos-express/marginal-eur/{chromosome}.txt"
  shell:
    """
    tar -Oxf {input} marginal_byChr_EURonly/EURonly_marginal_local_eQTL_meta_{wildcards.chromosome}.txt > \
    {output}
    """
    
rule filter_adipos_express_marginal_files:
  input: "data/adipos-express/marginal-eur/{chromosome}.txt"
  output: "data/adipos-express/processed-marginal-eur/{chromosome}.txt"
  shell: 
    """
    awk 'NR == 1 {{ print }} NR != 1 {{ if ($10 <= 5E-8) {{ print }} }} ' {input} > {output}
    """

rule process_adipos_express_marginal_data:
  input:
    adipos_express_marginal_paths = expand("data/adipos-express/processed-marginal-eur/{chromosome}.txt",
                                            chromosome = chromosomes),
    tss_data_path = "data/tss_data/hg19-tss-data.rds"                                      
  output:
    processed_data_path = "data/processed-data/adipos-express-marginal.rds"
  script: "code/process-data/process-adipos-express-marginal-data.R"
  
rule download_onek1k_data:
  output: "data/onek1k.tsv"
  shell:
    """
    wget https://onek1k.s3.ap-southeast-2.amazonaws.com/esnp/esnp_table.tsv.gz -O - | gunzip -c > {output}
    """
    
rule process_onek1k_data:
  input: 
    onek1k_path =  "data/onek1k.tsv",
    tss_data_path = "data/tss-data/hg19-tss-data.rds"
  output: processed_data_path = "data/processed-data/onek1k.rds"
  script: "code/process-data/process-onek1k-data.R"

rule download_metadata:
  output: 
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    pqtl_metadata_file = "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"
  shell:
    """
    wget -O {output.eqtl_metadata_file} https://zenodo.org/record/7808390/files/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz
    wget -O {output.pqtl_metadata_file} https://zenodo.org/record/7808390/files/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz
    """ 

rule download_tambets_etal_supp:
  output: "data/tambets-etal-supp.xlsx"
  shell:
    """
    wget -O {output} "https://www.biorxiv.org/content/biorxiv/early/2023/09/30/2023.09.29.560109/DC1/embed/media-1.xlsx?download=true"
    """

rule download_gnochhi_data: 
  output: "data/gnocchi-windows.bed" 
  shell: 
    """"
    wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06045-0/MediaObjects/41586_2023_6045_MOESM4_ESM.zip
    unzip 41586_2023_6045_MOESM4_ESM.zip 
    gunzip Supplementary_Data_2.bed.gz 
    mv Supplementary_Data_2.bed {output}
    rm 41586_2023_6045_MOESM4_ESM.zip 
    rm Supplementary_Data_1.tsv  Supplementary_Data_3.bed.gz Supplementary_Data_4.tsv
    rm Supplementary_Data_5.tsv Supplementary_Data_6_ESM.txt
    """

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

rule download_finngen_gwas_data: 
  output: "data/finngen/{gwas_id}.gz"
  localrule: True
  shell:
   """
   wget -O {output} https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_{wildcards.gwas_id}.gz
   """

rule process_finngen_gwas_data: 
  input: "data/finngen/{gwas_id}.gz"
  output: "data/finngen/{gwas_id}.gz.tbi"
  shell:
   """
   gunzip --force {input}
   bgzip data/finngen/{wildcards.gwas_id}
   tabix -s 1 -b 2 -e 2 -S 1 data/finngen/{wildcards.gwas_id}.gz
   """

rule download_finngen_manifest:
  output: "data/finngen/finngen-manifest.tsv"
  shell:
   """
   wget -O {output} https://storage.googleapis.com/finngen-public-data-r10/summary_stats/R10_manifest.tsv
   """

rule run_pqtl_eqtl_colocalisation:
  input: 
    eqtl_data_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz",
    eqtl_index_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz.tbi",
    permutation_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.permuted.tsv.gz",
    pqtl_data_file = "data/eqtl-catalogue/sumstats/QTD000584.cc.tsv.gz",
    pqtl_index_file = "data/eqtl-catalogue/sumstats/QTD000584.cc.tsv.gz.tbi",
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    pqtl_metadata_file = "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"
  output: 
    result_file = "output/data/pqtl-eqtl-coloc-{eqtl_id}-{chr}.rds"
  retries: 0
  resources: 
    mem_mb = lambda wildcards, attempt: 14000 * attempt,
    time_min = lambda wildcards, attempt: 120 * attempt
  script: "code/run-pqtl-eqtl-coloc-abf.R"

rule run_gwas_eqtl_colocalisation:
  input: 
    eqtl_data_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz",
    permutation_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.permuted.tsv.gz",
    eqtl_index_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz.tbi",
    gwas_data_file = "data/finngen/{gwas_id}.gz",
    gwas_index_file = "data/finngen/{gwas_id}.gz.tbi",
    eqtl_metadata_file = "data/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    manifest_file = "data/finngen/finngen-manifest.tsv"
  output: 
    coloc_results_file = "output/data/gwas-eqtl-coloc-{gwas_id}-{eqtl_id}-{chr}.rds",
    finemapping_results_file = "output/data/gwas-eqtl-finemapping-{gwas_id}-{eqtl_id}-{chr}.rds"
  retries: 0
  resources: 
    mem_mb = lambda wildcards, attempt: 14000 * attempt,
    time_min = lambda wildcards, attempt: 60 * attempt 
  script: "code/run-gwas-eqtl-coloc-abf.R"

rule make_ref_block:
  output: 
    leg_file = "data/{gene}.vcf.gz.impute.legend",
    haps_file = "data/{gene}.vcf.gz.impute.hap",
  script: "code/make-ref-block.sh"

rule run_simulations:
  input:
    leg_file = "data/{gene}.vcf.gz.impute.legend",
    haps_file = "data/{gene}.vcf.gz.impute.hap",
  output: sim_file = "output/data/sim-result-{gene}.rds"
  script: "code/run-simulation.R"
