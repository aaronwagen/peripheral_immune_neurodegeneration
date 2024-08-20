# Run coloc on Immune expression qtls and PD GWAS
# Aaron Wagen
# Jan 2023

# Run this script on a slurm cluster in parallel - one job for each eqtl (in this case for each chromosome per eqtl)
# Set up the script so the eqtl number is received from the submission command:
arguments <- commandArgs(trailingOnly = TRUE)
eqtl_number <- as.integer(arguments[1])
# eqtl_number <- 1

# Batch submission not working - it does not properly recognise the R project. This is a workaround, though it breaks the here::here usage.
setwd("/nemo/lab/gandhis/home/users/wagena/SNCA/")

# ml R/4.2.0-foss-2021b
# sbatch --time=3-00:00:00 --output logs/coloc_gwas_eqtl.log --open-mode=append --mem=128G --wrap="Rscript r_scripts/run_coloc_gwas_immune.R"

# Setup -------------------------------------------------------------------


library(here) # For project-specific paths
library(tidyverse) # For tidy manipulation of data
library(coloc)
library(colochelpR)
library(data.table)
library(plyranges)
library(tidytable) # For faster separation of columns of eqtl data
library(GenomicRanges) 
library(MafDb.1Kgenomes.phase3.GRCh38)
library(GenomicScores)
library(biomaRt)
library(qdapTools)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)


# Reading here
message("Reading here")
message(here::here())

message("Setting here")
#here::i_am("r_scripts/02i_run_coloc_gwas_immune_aquino_rizig.R")

here::here()

# Assign arguments
args <- list(
  gwas_name = "rizig2020",
  gwas_dir = here::here("raw_data/GWAS/coloc_munged_data"),
  munged_gwas = "/nemo/lab/gandhis/home/users/wagena/references/GWAS/coloc_munged_data/rizig2023_gwas_tidy_varbeta.txt",
  PD_gwas_details = list(path = "/nemo/lab/gandhis/home/users/wagena/references/GWAS/Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38.txt",
                     n_cases = 1488,
                     n_controls = 196430,
                     n_total = 197918,
                     prop_cases = 0.007518265,
                     cc_or_quant = "cc"),
  eqtl_dir = here::here("raw_data/GWAS/aquino_2023_covid_eqtl"),
  eqtl_sample_number = 222, # There are 222 individuals included in the aquino dataset
  output_path = here::here("processed_data",
                           "rizig_aquino",
                           "coloc"),
  ensembl_gene_ids_overlapping_1Mb_window_hit_path =  here::here("raw_data/GWAS/coloc_munged_data/rizig_ensembl_gene_ids_overlapping_1Mb_window_hit.RDS"),
  gtf_38.84 = "/nemo/lab/gandhis/home/users/wagena/references/GRCh38/annotation/Homo_sapiens.GRCh38.84.gtf.gz"
)


gwas_results_path = here::here("raw_data/GWAS/coloc_munged_data")

# Note: here using the 1000 genomes project east asian population as it gave more results for MAF than other databases (eg gnomad 2. gnomad3 MafH5.gnomAD.v3.1.2.GRCh38 did not have an east asian population)
mafdb1K <- MafDb.1Kgenomes.phase3.GRCh38

# Load pretidied GWAS data ---------------------------------------------------------------


# GWAS

print("Loading GWAS")
gwas_tidy <- read_delim(args$munged_gwas)


# The GWAS has been previously parsed and tidied as per 01_munge_GWASs_for_coloc.R. 
# For aquino we are using the rs number as the SNP identifier


# Find all genes within +/-1Mb of significant GWAS hits -------------------------

print("Extracting genes of interest")

# Extract all genes within +/- 1Mb of significant RBD hits witht the following code. I have pre-saved the output as the ensembl portal often fails.
# ensembl_gene_names_to_test  <- gwas_tidy %>%
#   dplyr::mutate(CHR = as.integer(CHR)) %>%
#   colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
#                                                   CHR_column = "CHR",
#                                                   BP_column = "BP",
#                                                   mart = 38)

ensembl_gene_names_to_test  <- readRDS(args$ensembl_gene_ids_overlapping_1Mb_window_hit_path)  


# Load aquino eqtl snp information ------------------------------------------

# In this dataset:
## eQTLs were run in 4 different conditions:
## 1) **celltype_condition___CellPropLineage_SVs_RUN1**
##   - Expression QTL mapping at cell type resolution (ie, resolution level 2), adjusted on proportions within the lineage and SVs
## 2) **celltype_condition_logFC__CellPropLineage_SVs_RUN1**
##   - Response QTL mapping at cell type resolution (ie, resolution level 2), adjusted on proportions within the lineage and SVs
## 3) **lineage_condition___CellPropLineage_SVs_RUN1**
##   - Expression QTL mapping at lineage resolution (ie resolution level 1) adjusted on proportions within the lineage and SVs
## 4) **lineage_condition_logFC__CellPropLineage_SVs_RUN1**
##   - Response QTL mapping at lineage resolution (ie resolution level 1), adjusted on proportions within the lineage and SVs

## A1 allele is the allele in the reference genome - this is coloc Al2
## A2 allele is the effect allele - this is coloc Al1
## The p value is a raw p value. We will manually FDR these values
## There are blood donations from 982 participants. The various blood cells were isolated from these samples so each eqtl has 416 samples.

print("Load eqtl")

## Select eqtl
# 
# eqtl_list <- list.files(args$eqtl_dir, 
#                         pattern = "^eQTL_FineMapped_.*\\.txt\\.gz$", 
#                         full.names = TRUE, 
#                         recursive = TRUE)
# 
# 
# eqtl_df <- tibble(eqtl_path = eqtl_list) %>% 
#   dplyr::mutate(eqtl_id = str_remove(eqtl_path, "^/nemo/lab/gandhis/home/users/wagena/SNCA/raw_data/GWAS/aquino_2023_covid_eqtl/")) %>%
#   separate(eqtl_id, into = c("cell_lineage_condition", "celltype_virus", "chromosome_file"), sep = "/") %>%
#   separate(celltype_virus, into = c("celltype", "virus"), sep = "__") %>%
#   dplyr::mutate(
#     eqtl_chromosome = str_extract(chromosome_file, "chr[0-9]+"),
#     eqtl_resolution_level = case_when(
#       str_detect(cell_lineage_condition, "lineage_condition") ~ "celltype",
#       str_detect(cell_lineage_condition, "celltype_condition") ~ "subtype"),
#     expression_response = case_when(
#       str_detect(cell_lineage_condition, "logFC") ~ "response",
#       TRUE ~ "expression"))
# 
# 
# cell_names <- tibble(celltype = eqtl_df$celltype %>% unique() %>% sort(),
#                      cell_name = c("B-cell",
#                                    "B-memory-k",
#                                    "B-memory-l",
#                                    "B-naive-k",
#                                    "B-naive-l",
#                                    "Dendritic-classic",
#                                    "NK-innate_lymphoid",
#                                    "T-mucosal-assoc-invariant",
#                                    "Monocyte",
#                                    "Mono-cd14",
#                                    "Mono-cd14-infected",
#                                    "Mono-cd16",
#                                    "NK",
#                                    "NK-CD56-bright",
#                                    "NK-CD56-dim",
#                                    "NK-memory-like",
#                                    "Dendritic-plasmacytoid",
#                                    "Plasmablast",
#                                    "T-CD4",
#                                    "T-CD4-effector",
#                                    "T-CD4-naive",
#                                    "T-CD8",
#                                    "T-CD8-central-effector-mem",
#                                    "T-CD8-effector-mem-reexpress-cd45RA",
#                                    "T-CD8-naive",
#                                    "T-gammadelta",
#                                    "T-reg")
# )
# 
# eqtl_df <- left_join(eqtl_df,
#                      cell_names,
#                      by = "celltype") %>% 
#   dplyr::mutate(eqtl_name = str_c(expression_response, "_", eqtl_resolution_level, "_", cell_name, "_", virus )) %>% 
#   dplyr::select(eqtl_name, eqtl_path, celltype, virus, expression_response, cell_name, eqtl_chromosome)
# 
# 
# write_csv(eqtl_df,
#           file = str_c(args$eqtl_dir, "/aquino_summary_table.csv"))

eqtl_df <- read_csv(file = str_c(args$eqtl_dir, "/aquino_summary_table.csv"))


# Coloc loop --------------------------------------------------------------


# eqtl_number <- 1
print(eqtl_number)

eqtl_name <- eqtl_df$eqtl_name[eqtl_number] 
eqtl_path <- eqtl_df$eqtl_path[eqtl_number]
eqtl_sample_no <- args$eqtl_sample_number


print(str_c("eQTL is:   ", eqtl_name))
print(str_c("eqtl Path is:   ", eqtl_path))
print(str_c("Number of samples per tissue is: ", eqtl_sample_no))


## Load eqtl
eqtl <- fread(eqtl_path) # Load eqtl


eqtl <- eqtl %>% # Perform minimal munging on overall qtl here to allow for munging at a smaller scale within loops below.
  dplyr::mutate(FDR = p.adjust(pvalue, method = "fdr")) %>% 
  dplyr::select(SNP = rsID, CHROM, POS_B38, gene, Al2 = REF, Al1 = ALT, celltype, condition, beta = slope, se = slope_se, p.value = FDR) 


# Find overlapping genes from qtl and gwas hits to save processing time in the next loop
eqtl_genes <- unique(eqtl$gene)
cis_genes_matched_gwas_eqtl <- intersect(eqtl_genes, ensembl_gene_names_to_test)

## Create results directory
results_path_GWAS_region_eqtl <- colochelpR::make_results_dir(results_path = args$output_path, folder_name = eqtl_name)

## Within results directory, create a folder for "liberal" and "robust" coloc p12 prior
results_path_priors <- setNames(vector(mode = "list", length = 2),
                                c("liberal", "robust"))

results_path_priors$liberal <- make_results_dir(results_path = results_path_GWAS_region_eqtl, folder_name = "liberal")
results_path_priors$robust <- make_results_dir(results_path = results_path_GWAS_region_eqtl, folder_name = "robust")


## Set p value cuttofs
p12 <- setNames(c(1e-05, 5e-06),
                c("liberal", "robust"))



for(j in seq_along(cis_genes_matched_gwas_eqtl)){
  
  # j = 1
  ensembl_gene_id_to_search <-  cis_genes_matched_gwas_eqtl[j]
  
  
  print(str_c(Sys.time(), " - foo_gwas - ", eqtl_name, "_eqtl - ", ensembl_gene_id_to_search))
  
  # Filter eqtls for snps matching a single editing site
  eqtl_tidy_gene_filtered <- 
    eqtl %>% 
    dplyr::filter(gene == ensembl_gene_id_to_search)    %>% 
    dplyr::filter(nchar(Al1) == 1,  # Remove the SNPs in the eqtl that are insertions/deletions (ie not true SNPs)
                  nchar(Al2) == 1) %>% 
    dplyr::mutate(eqtl_dataset = eqtl_name,
                  N = eqtl_sample_no) %>% 
    as_granges(seqnames = CHROM,    # Make genomic range to get mafs.
               start = POS_B38,
               end = POS_B38) %>% 
    gscores(x = mafdb1K, ranges = ., pop = "AF") %>% # This is a mixed population so use "AF" for overall allele frequence
    as_tibble() %>% 
    dplyr::rename(maf = AF) %>% 
    dplyr::select(eqtl_dataset,
                  gene,
                  SNP,
                  beta,
                  se,
                  p.value,
                  Al1,
                  Al2,
                  maf,
                  N) %>% 
    dplyr::filter(!is.na(maf), #Remove NAs and 0 from maf and beta, 
                  maf != 0,
                  beta != 0) %>%
    dplyr::mutate(maf = dplyr::case_when(maf > 0.5 ~ 1-maf, # flip mafs so they are all <0.5
                                         TRUE ~ maf)) %>% 
    colochelpR::get_varbeta() %>% 
    colochelpR::check_coloc_data_format(.,
                                        beta_or_pval = "beta",
                                        check_maf = T) %>%
    dplyr::group_by(SNP) %>% # Remove duplicates choosing the one with the lowest p value
    arrange(p.value) %>%  # Arrange rows within each SNP group by p-value
    distinct(SNP, .keep_all = TRUE)
  
  
  if (nrow(eqtl_tidy_gene_filtered) == 0) {
    print(str_c("No QTLs overlapping: ", ensembl_gene_id_to_search))
    next
  }
  
  # Run coloc
  p12 <- setNames(c(1e-05, 5e-06),
                  c("liberal", "robust"))
  
  for(k in 1:length(p12)){
    
    print(str_c("Results for '", names(p12[k]), "' p12 prior;  p12 = ", p12[k]))
    
    coloc_results_annotated <-
      colochelpR::get_coloc_results(df1 = gwas_tidy, df2 = eqtl_tidy_gene_filtered, 
                                    # Harmonise set to true as is will flip b-values to allow exploration of directionality of results
                                    harmonise = T, 
                                    df1_type = "cc", df2_type = "quant", 
                                    df1_beta_or_pval = "beta", df2_beta_or_pval = "beta",
                                    df_1_propor_cases = args$PD_gwas_details$prop_cases,
                                   # df1_N = max(gwas_tidy$N), # df_N not required for case-control datasets
                                    df2_N = max(eqtl_tidy_gene_filtered$N),
                                    annotate_signif_SNP_df1_df2 = T, 
                                    key_cols = c("GWAS_1", "eqtl_dataset_2", "gene_2"), 
                                    df_1_name = "GWAS", df_2_name = "eqtl", 
                                    df1_path = args$gwas_dir, df2_path = args$eqtl_dir,
                                    p1 = 1e-04, p2 = 1e-04, p12 = as.numeric(p12[k]))
    
    colochelpR::save_coloc_results(coloc_results_annotated, results_dir_path = results_path_priors[[names(p12[k])]])
  }
  
}



print("Done")
