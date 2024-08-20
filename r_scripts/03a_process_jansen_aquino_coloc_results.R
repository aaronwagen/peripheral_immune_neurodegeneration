#
# Title: "Analyse colocalisation of peripheral immune-cell genes and PD GWAS risk in the foo_aquino gwas"
# Author: "Aaron Wagen"
# This script parses the output folders of the run_coloc_gwas_immune.R script and generates a summary table of results

# Run this script on a slurm cluster with:
# ml R/4.2.0-foss-2021b
# sbatch --time=8:00:00 --output logs/coloc_jansen_aquino_processing.log --open-mode=truncate --mem=64G --wrap="Rscript r_scripts/03k_process_jansen_aquino_coloc_results.R"


## Load libraries
library(here) # For project-specific paths
library(tidyverse) # For tidy manipulation of data
library(coloc)
library(colochelpR)
library(data.table)
library(plyranges)
library(tidytable) # For faster separation of columns of eqtl data
library(GenomicRanges) 
library(biomaRt)
library(qdapTools)
library(LDlinkR)
library(pheatmap)
library(DBI)
library(patchwork)
library(viridis)
library(sciRmdTheme)
library(paletteer)

# Set path
here::i_am("r_scripts/03k_process_jansen_aquino_coloc_results.R")

args <- list(
  # munged_gwas = "/nemo/lab/gandhis/home/users/wagena/references/GWAS/coloc_munged_data/PD2019_ex23andMe_tidy_varbeta.txt",
  # eqtl_dir = here::here("raw_data", 
  #                       "ebi_eqtl_catalogue"),
  # eqtl_sample_info_file = here::here("raw_data", 
  #                                    "ebi_eqtl_catalogue",
  #                                    "ebi_eqtl_metadata.tsv"),
  coloc_path = here::here("processed_data",
                           "jansen_aquino",
                           "coloc"),
  coloc_results_summary_file = here::here("processed_data",
                                          "jansen_aquino",
                                          "coloc",
                                          "coloc_summary.Rds"),
  gencode_v41_gtf = "/camp/home/wagena/references/GRCh38/annotation/gencode.v41.annotation.gtf"
  
)


# Load biomart function with mirror

biomart_df <- function (dataframe, 
          columnToFilter, 
          mart = 38, 
          attributes, 
          filter,
          mirror = "asia") {
  if (mart != 38 && mart != 37) 
    stop("Mart must be 38 or 37...")
  if (mart == 38) {
    ensembl_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                     dataset = "hsapiens_gene_ensembl",
                                     host = str_c("https://", mirror, ".ensembl.org"))
  }
  else if (mart == 37) {
    ensembl_mart <- biomaRt::useMart(host = "grch37.ensembl.org", 
                                     biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  }
  genes <- dataframe %>% .[[columnToFilter]] %>% unique()
  print(stringr::str_c("Number of unique genes to search: ", 
                       length(genes)))
  biomart_query <- biomaRt::getBM(attributes = attributes, 
                                  filters = filter, values = genes, mart = ensembl_mart)
  print(stringr::str_c("Number of matches found:", nrow(biomart_query)))
  join_vector <- filter
  names(join_vector) <- columnToFilter
  merged <- dplyr::inner_join(dataframe, biomart_query, by = join_vector)
  return(merged)
}


# Methods


# Process Coloc results


# Create table with file paths, datasets, prior type and editing site

results_dir <- tibble(file_path =  list.files(args$coloc_path, recursive = T, full.names = T, pattern = ".rda")) %>% 
  dplyr::mutate(dir = file_path %>% 
                  str_replace(., "/[^/]*$", "") %>% 
                        str_replace(., "/[^/]*$", ""),
                file_name = basename(file_path) %>% 
                         str_replace(., ".rda", ""),
                dataset = basename(dir),
                prior_type = file_path %>% 
                  str_replace(., "/[^/]*$", "") %>% 
                  basename(),
                gene = str_replace(file_name, ".+_(.*)", "\\1") # Select out everything after the last hyphen  - the gene name
                ) 


# Split into liberal/robust results
dataset_dir <- setNames(results_dir %>% dplyr::group_split(prior_type),
                        c(results_dir %>% .[["prior_type"]] %>% unique() %>% sort()))
# Priors df
priors_df <- tibble(prior_type = c("liberal", "robust"),
                    p12 = c(1e-05, 5e-06))
all_results <- setNames(vector(mode = "list", length = length(dataset_dir)),
                        names(dataset_dir))


# List of datasets analysed
dataset_list <- unique(results_dir$dataset)



for(i in 1:length(dataset_dir)){
  
  priors_df_filtered <- 
    priors_df %>% 
    dplyr::filter(prior_type == names(dataset_dir[i]))
  
  dataset_file_paths <- dataset_dir[[i]]
  
  dataset_names <- dataset_file_paths$dataset %>% unique()
  
  results_list <- vector(mode = "list", length = length(dataset_names))
  
  for(j in 1:length(dataset_names)){
    
    dataset <- dataset_names[j]
    
    print(dataset)
    
    dir_to_load <- 
      dataset_file_paths %>% 
      dplyr::filter(dataset == dataset_names[j]) %>% 
      .[["dir"]] %>% 
      unique()
    
    dir_to_load <- str_c(dir_to_load, "/", priors_df_filtered$prior_type)
    
    for(k in 1:length(dir_to_load)){
      
      print(dir_to_load[k])

      if(dataset %in% dataset_list){
        
        results <- 
          colochelpR::merge_coloc_summaries(dir_to_load[k], add_signif_SNP = F, recursive = T, pattern = ".rda") %>% 
          dplyr::select(GWAS_1, gene_2, everything(), -eqtl_dataset_2)
        
      }
      
      results <- 
        results %>% 
        dplyr::mutate(dataset = dataset,
                      p12 = priors_df_filtered$p12) %>% 
        dplyr::select(GWAS_1, dataset, gene_2, everything())
      
      if(k == 1){
        
        results_list[[j]] <- results 
        
      } else{
        
        results_list[[j]] <- 
          results_list[[j]] %>% 
          dplyr::bind_rows(results)
        
      }
      
    }
    
  }
  
  all_results[[i]] <- results_list
  
}


# Biomart is rubbish so use manual gtf to find gene names. Matching gencode v27 which is used in the eqtl.
gtf <- rtracklayer::readGFF(args$gencode_v41_gtf) %>% 
  as_tibble() %>% 
  dplyr::select(gene_id, gene_name) %>% 
  dplyr::distinct(gene_name, gene_id) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

results <- 
  all_results %>% 
  lapply(., function(x){
    
    x %>% 
      qdapTools::list_df2df() %>% 
      left_join(.,
                gtf,
                by = c("gene_2" = "gene_id")) %>% 
      dplyr::select(GWAS_1, dataset, gene_2, everything(), -X1)
    
  }) 


saveRDS(results, file = args$coloc_results_summary_file)
  

