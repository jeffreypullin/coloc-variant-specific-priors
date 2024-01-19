configfile: "config.yaml"

chromosomes = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]

# .cc.tsv.gz
eqtl_data_ids_subset = ["QTD000110", "QTD000373", "QTD000356"]

rule all: 
  input: 
    "data/processed-data/eqtl-catalogue.rds",
    "data/processed-data/gtex.rds",
    "data/processed-data/adipos-express-marginal.rds",
    "data/processed-data/eqtlgen.rds",
    "data/processed-data/onek1k.rds",
    "data/fauman-hyde/eqtlgen.txt",
    "data/tss-data/hg19-tss-data.rds",
    expand("data/adipos-express/processed-ab1-eur/{chromosome}.txt", chromosome = chromosomes),
    expand("output/data/pqtl-eqtl-coloc-{eqtl_id}-{chr}.tsv", 
           chr = [x for x in range(1, 23)],
           eqtl_id = eqtl_data_ids_subset)

rule download_tss_data: 
  output: hg19_tss_data_path = "data/tss-data/hg19-tss-data.rds",
          hg38_tss_data_path = "data/tss-data/hg38-tss-data.rds"
  script: "code/download-tss-data.R"
    
rule download_fauman_hyde_eqtlgen_data: 
  output: "data/fauman-hyde/eqtlgen.txt"
  shell:
    """
    wget -O {output} https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-022-04706-x/MediaObjects/12859_2022_4706_MOESM3_ESM.txt
    """

rule download_eqtlgen_data:
  output: "data/eqtlgen.txt"
  shell:
    """
    wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
    gunzip -d 2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
    mv 2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt {output}
    """
    
rule process_eqtlgen_data:
  input: eqtlgen_path =  "data/eqtlgen.txt"
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
    tabix_file =  "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz.tbi"
  script: "code/download-eqtl-catalogue.R"
  
rule filter_eqtl_catalogue_files:
  input: "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz"
  output: "data/eqtl-catalogue/processed-sumstats/{dataset_id}.cc.tsv"
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

rule run_pqtl_eqtl_colocalisation:
  input: 
    eqtl_data_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz",
    eqtl_index_file = "data/eqtl-catalogue/sumstats/{eqtl_id}.cc.tsv.gz.tbi",
    pqtl_data_file = "data/QTD000584.cc.tsv.gz",
    pqtl_index_file = "data/QTD000584.cc.tsv.gz.tbi",
    eqtl_metadata_file = "data/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
    pqtl_metadata_file = "data/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"
  output: 
    result_file = "output/data/pqtl-eqtl-coloc-{eqtl_id}-{chr}.tsv"
  script: "code/run-coloc-abf.R"
