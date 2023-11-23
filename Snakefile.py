configfile: "config.yaml"

chromosomes = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]

rule all: 
  input: 
    "data/plot-data/eqtl-catalogue.rds",
    "data/plot-data/gtex.rds",
    "data/plot-data/adipos-express-marginal.rds",
    "data/plot-data/eqtlgen.rds",
    "data/plot-data/onek1k.rds",
    "data/fauman-hyde/eqtlgen.txt",
    expand("data/adipos-express/processed-ab1-eur/{chromosome}.txt", chromosome = chromosomes),

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
    
rule extract_plot_data_eqtlgen:
  input: eqtlgen_path =  "data/eqtlgen.txt"
  output: plot_data_path = "data/plot-data/eqtlgen.rds"
  script: "code/extract-plot-data-eqtlgen.R"
    
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
    
rule extract_plot_data_gtex: 
  input: gtex_paths = expand("data/gtex-v8/{tissue}.v8.signif_variant_gene_pairs.txt", tissue = config["gtex_tissues"]), 
  output: plot_data_path = "data/plot-data/gtex.rds"
  script: "code/extract-plot-data-gtex.R"

rule download_etl_catalogue_metadata: 
  output: "data/eqtl-catalogue/eqtl-catalogue-metadata.tsv"
  shell:
    """
    wget -O {output} https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv
    """

rule download_eqtl_catalogue_files: 
  input: metadata_file = "data/eqtl-catalogue/eqtl-catalogue-metadata.tsv"
  output: dataset_file = "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz"
  script: "code/download-eqtl-catalogue.R"
  
rule process_eqtl_catalogue_files:
  input: "data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz"
  output: "data/eqtl-catalogue/processed-sumstats/{dataset_id}.cc.tsv"
  shell:
    """
    zcat {input} | awk 'NR == 1 {{ print }} NR != 1 {{ if ($9 <= 5E-8) {{ print }} }} ' > {output}
    """

rule extract_plot_data_eqtl_catalogue:
  input:
    eqtl_catalogue_metadata = "data/eqtl-catalogue/eqtl-catalogue-metadata.tsv",
    eqtl_catalogue_paths = expand("data/eqtl-catalogue/processed-sumstats/{dataset_id}.cc.tsv",
                                  dataset_id = config["eqtl_catalogue_dataset_ids"]),
  output:
    plot_data_path = "data/plot-data/eqtl-catalogue.rds"
  script: "code/extract-plot-data-eqtl-catalogue.R"
  
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
    
rule process_adipos_express_marginal_files:
  input: "data/adipos-express/marginal-eur/{chromosome}.txt"
  output: "data/adipos-express/processed-marginal-eur/{chromosome}.txt"
  shell: 
    """
    awk 'NR == 1 {{ print }} NR != 1 {{ if ($10 <= 5E-8) {{ print }} }} ' {input} > {output}
    """

rule extract_plot_data_adipos_express_marginal:
  input:
    adipos_express_marginal_paths = expand("data/adipos-express/processed-marginal-eur/{chromosome}.txt",
                                            chromosome = chromosomes)
  output:
    plot_data_path = "data/plot-data/adipos-express-marginal.rds"
  script: "code/extract-plot-data-adipos-express-marginal.R"
  
rule download_onek1k_data:
  output: "data/onek1k.tsv"
  shell:
    """
    wget https://onek1k.s3.ap-southeast-2.amazonaws.com/esnp/esnp_table.tsv.gz -O - | gunzip -c > {output}
    """
    
rule extract_plot_data_onek1k:
  input: onek1k_path =  "data/onek1k.tsv"
  output: plot_data_path = "data/plot-data/onek1k.rds"
  script: "code/extract-plot-data-onek1k.R"
