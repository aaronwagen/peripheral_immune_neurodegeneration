# Run coloc on Immune expression qtls and PD GWAS
# Aaron Wagen
# Jan 2023

# Run this script on a slurm cluster with:
# ml R/4.2.0-foss-2021b
# sbatch --time=3-00:00:00 --output logs/coloc_rizig2023_nedelecquach.log --open-mode=truncate --mem=128G --wrap="Rscript r_scripts/02c_run_coloc_gwas_immune_african.R"

# Setup -------------------------------------------------------------------


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
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)
library(MafDb.1Kgenomes.phase3.GRCh38)
library(GenomicScores)
library(parallel)
library(foreach)



here::i_am("r_scripts/02c_run_coloc_gwas_immune_african.R")



# GWAS datasets:

## 3) Rizig Bandres-Ciga et al, medRxiv, 2023, Genome-wide Association Identifies Novel Etiological Insights Associated with Parkinson’s Disease in African and African Admixed Populations
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10312852/
## Columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, variant_id, ref_allele, direction, HetISq, HetChiSq, HetDf, HetPVal, rsid
## This includes 1488 PD cases, and 196,430 controls
## File supplied by authors via C Blauwendraat.


# eQTL datasets:

## 


# Assign arguments
args <- list(
  munged_gwas = here::here("raw_data",
                           "rizig2023_gwas_tidy_varbeta.txt"),
  eqtl_dir = here::here("raw_data", 
                         "ebi_eqtl_catalogue"),
  eqtl_sample_info_file = here::here("raw_data", 
                                  "ebi_eqtl_catalogue",
                                  "ebi_eqtl_metadata.tsv"),
  output_path = here::here("processed_data",
                           "rizig_2023",
                           "coloc"),
  gwas_details = list(n_cases = 1488,
                      n_controls = 196430,
                      n_total = 197918,
                      prop_cases = 0.007518265,
                      cc_or_quant = "cc")
)




# Load pretidied GWAS data ---------------------------------------------------------------

# GWAS


dir.create(args$output_path)

print("Loading GWAS")
gwas_tidy <- read_delim(args$munged_gwas)



# Find all genes within +/-1Mb of significant GWAS hits -------------------------

print("Extracting genes of interest")

# Use the following code to extract all genes within +/- 1Mb of significant GWAS hits with the following code.
ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas_tidy %>%
  colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                                  CHR_column = "CHR",
                                                  BP_column = "BP",
                                                  mart = 38)
saveRDS(ensembl_gene_ids_overlapping_1Mb_window_hit,
        file = str_c(args$output_path, "/ensembl_gene_ids_overlapping_1Mb_window_hit.RDS")) # Save this to a file to prevent errors when trying to connect to biomart
ensembl_gene_names_to_test <- readRDS(str_c(args$output_path, "/ensembl_gene_ids_overlapping_1Mb_window_hit.RDS"))


# Use the following code to look at specific genes
# ensembl_gene_names_to_test <- c("LRRK2",   "BST1",    "MCCC1",  "VPS13C",  "DGKQ", "SNCA", "TMEM175") # These have already been processed
# ensembl_gene_names_to_test <- c("MMRN1", "SNCA-AS1") # Additional genes to run


# Load oto et al eqtl snp information ------------------------------------------

# In this dataset:
## Ref allele is the allele in the reference genome - this is coloc Al2
## Alt allele is the effect allele - this is coloc Al1
## There are blood donations from 416 participants. The various blood cells were isolated from these samples so each eqtl has 416 samples.

print("Load eqtl")

## Select eqtl

eqtl_list <- list.files(args$eqtl_dir, pattern = "*.tsv.gz")  # Load eqtls

eqtl_metadata <- read_delim(args$eqtl_sample_info_file) %>% 
  dplyr::mutate(eqtl_name = str_c(study_label, "_", sample_group, "_", dataset_id),
                eqtl_file = str_c(eqtl_name, ".all.tsv.gz")) %>% 
  dplyr::filter(eqtl_file %in% eqtl_list)


# Coloc loop --------------------------------------------------------------


cl <- parallel::makeCluster(detectCores()) # Detect number of cores by which to parallelise
doParallel::registerDoParallel(cl) # Set up number of clusters

for(eqtl_number in 1:4) {
  
  print(eqtl_number)
  
  eqtl_sample <- eqtl_metadata$eqtl_file[eqtl_number]
  eqtl_name <- eqtl_metadata$eqtl_name[eqtl_number] %>% 
    str_replace("_QTD.*", "")
  eqtl_path <- str_c(args$eqtl_dir, "/", eqtl_sample)
  eqtl_sample_no <- eqtl_metadata$sample_size[eqtl_number]
  # eqtl_sample_no <- eqtl_sample_numbers %>%  # Retreive number of samples in eqtl
  #   dplyr::filter(tissue == eqtl_name) %>%
  #   dplyr::select(sample_no) %>%
  #   unlist()
  
  print(str_c("Tissue is:   ", eqtl_sample))
  print(str_c("Tissue Name is:   ", eqtl_name))
  print(str_c("Tissue Path is:   ", eqtl_path))
  print(str_c("Number of samples per tissue is: ", eqtl_sample_no))
  
  
  ## Load eqtl
  eqtl <- fread(eqtl_path) # Load eqtl
  
  eqtl <- eqtl %>% # Perform minimal munging on overall qtl here to allow for munging at a smaller scale within loops below.
    dplyr::select(-c(type, r2, ac, an, median_tpm, gene_id, molecular_trait_object_id))
  
  # Find overlapping genes from qtl and gwas hits to save processing time in the next loop
  eqtl_genes <- unique(eqtl$molecular_trait_id)
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
    
    ensembl_gene_id_to_search <-  cis_genes_matched_gwas_eqtl[j]
    
    
    print(str_c(Sys.time(), " - PD_gwas - ", eqtl_name, "_eqtl - ", ensembl_gene_id_to_search))
    
    # Filter eqtls for snps matching a single editing site
    eqtl_tidy_gene_filtered <- 
      eqtl %>% 
      dplyr::filter(molecular_trait_id == ensembl_gene_id_to_search)    %>% 
      dplyr::filter(nchar(ref) == 1,  # Remove the SNPs in the eqtl that are insertions/deletions (ie not true SNPs)
                    nchar(alt) == 1) %>% 
      dplyr::mutate(eqtl_dataset = eqtl_name,
                    N = eqtl_sample_no) %>% 
      dplyr::rename(gene = molecular_trait_id,
                    SNP = rsid,
                    p.value = pvalue,
                    Al1 = alt,
                    Al2 = ref) %>% 
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
                                      df1_beta_or_pval = "beta", df2_beta_or_pval = "pval",
                                      df_1_propor_cases = args$gwas_details$prop_cases,
                                      #df1_N = args$gwas_details$n_total # df_N not required for case-control datasets
                                      df2_N = max(eqtl_tidy_gene_filtered$N),
                                      annotate_signif_SNP_df1_df2 = T, 
                                      key_cols = c("GWAS_1", "eqtl_dataset_2", "gene_2"), 
                                      df_1_name = "GWAS", df_2_name = "eqtl", 
                                      df1_path = args$gwas_dir, df2_path = args$eqtl_dir,
                                      p1 = 1e-04, p2 = 1e-04, p12 = as.numeric(p12[k]))
      
      colochelpR::save_coloc_results(coloc_results_annotated, results_dir_path = results_path_priors[[names(p12[k])]])
      

    }
    
  }
  
  rm(eqtl) # To save memory before loading the next eqtl
}


stopImplicitCluster()  # Stop clusters

print("Done")
