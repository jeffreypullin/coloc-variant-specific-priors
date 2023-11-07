configfile: "config.yaml"

chromosomes = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]
endothelial_celltypes = ["Lymphatic", "Systemicvenous", "aCap", "arteriole", "gCap", "venule"]

rule all: 
  input: 
    "data/fauman-hyde/eqtlgen.txt",
    "data/eqtlgen.txt",  
    expand("data/gtex-v8/{tissue}.v8.signif_variant_gene_pairs.txt", tissue = config["gtex_tissues"]), 
    expand("data/eqtl-catalogue/sumstats/{dataset_id}.cc.tsv.gz", dataset_id = config["eqtl_catalogue_dataset_ids"]),
    expand("data/adipos-express/marginal-eur/{chromosome}.txt", chromosome = chromosomes),
    expand("data/adipos-express/ab1-eur/{chromosome}.txt", chromosome = chromosomes),
    expand("data/lung-single-cell/endothelial/{celltype}.txt", celltype = endothelial_celltypes)

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

rule download_etl_catalogue_metadata: 
  output: "data/eqtl-catalogue/" + "eqtl-catalogue-metadata.tsv"
  shell:
    """
    wget -O {output} https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv
    """

rule download_eqtl_catalogue_files: 
  input: 
    metadata_file = "data/eqtl-catalogue/" + "eqtl-catalogue-metadata.tsv"
  output: 
    dataset_file = "data/eqtl-catalogue/sumstats/" + "{dataset_id}.cc.tsv.gz"
  script: "code/download-eqtl-catalogue.R"
  
rule download_adipos_express_tars: 
  output: 
    marginal = temporary("data/adipos-express-marginal-eur-by-chr.tar"),
    ab1 = temporary("data/adipos-express-ab1-eur-by-chr.tar")
  shell:
    """
    wget ftp://mohlkeanon:anon@rc-ns-ftp.its.unc.edu/marginal_byChr_EURonly.tar.gz -O - | gunzip -d > {output.marginal}
    wget ftp://mohlkeanon:anon@rc-ns-ftp.its.unc.edu/AB1_byChr_EURonly.tar.gz -O - | gunzip -d > {output.ab1}
    """

rule extract_adipos_express_marginal_files:
  input: "data/adipos-express-marginal-eur-by-chr.tar",
  output: "data/adipos-express/marginal-eur/{chromosome}.txt"
  shell:
    """
    tar -Oxf {input} marginal_byChr_EURonly/EURonly_marginal_local_eQTL_meta_{wildcards.chromosome}.txt > \
    {output}
    """
    
rule extract_adipos_express_ab1_files:
  input: "data/adipos-express-ab1-eur-by-chr.tar",
  output: "data/adipos-express/ab1-eur/{chromosome}.txt"
  shell:
    """
    tar -Oxf {input} AB1_byChr_EURonly/EURonly_AB1_local_eQTL_meta_{wildcards.chromosome}.txt > \
    {output}
    """

rule download_lung_endothelial_tar:
  output: temporary("data/endothelial.tar")
  shell:
    """
    wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227136/suppl/GSE227136%5Flimix%5Fres%5Fendothelial.tar.gz" -O - | gunzip -d > {output}
    """
    
rule extract_lung_endothelial_files:
  input: "data/endothelial.tar"
  output: "data/lung-single-cell/endothelial/{celltype}.txt"
  shell:
    """
    tar -Oxf {input} endothelial_{wildcards.celltype}_qtl_results_all.txt > {output}
    """