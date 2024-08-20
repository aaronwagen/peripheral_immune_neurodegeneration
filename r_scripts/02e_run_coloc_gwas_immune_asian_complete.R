# Run coloc on Immune expression qtls and PD GWAS
# Aaron Wagen
# Jan 2023

# Run this script on a slurm cluster with:
# ml R/4.2.0-foss-2021b
# sbatch --time=2-00:00:00 --output logs/coloc_foo2020_ota2021_complete.log --open-mode=truncate --mem=128G --wrap="Rscript r_scripts//02b_run_coloc_gwas_immune_asian_complete.R"

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



here::i_am("r_scripts/02b_run_coloc_gwas_immune_asian.R")



# GWAS datasets:

## 1) Sakaue and Kanai - A cross population atlas of genetic associations for 220 human phenotypes. 
## https://www.nature.com/articles/s41588-021-00931-x#Sec27
## This includes 340 PD cases, and 175788 controls.
## Data downloaded from: https://humandbs.biosciencedbc.jp/en/hum0197-v3-220
## Dictionary file at: https://humandbs.biosciencedbc.jp/files/hum0197.org/Dictfile_BBJ.html


## 2) Foo, Chew, Chung et al - Identification of Risk Loci for Parkinson Disease in Asians and Comparison of Risk Between Asians and Europeans
## https://jamanetwork.com/journals/jamaneurology/fullarticle/2764340
## columns: CHR, BP, SNP (rs number:location:A1:A2), A1, A2, BETA, P SE, MarkerName (rs number)
## This includes 6724 PD cases, and 24581 controls, 5,843,213 snps
## File supplied by authors via C Blauwendraat.
## This has previously been munged as per 01_munge_GWASs_for_coloc.R




# Assign arguments
args <- list(
  munged_gwas = here::here("raw_data",
                           "foo2020_gwas_tidy_varbeta.txt"),
  gwas_details = list(n_cases = 6724,
                      n_controls = 24581,
                      n_total = 31305,
                      prop_cases = 0.214790,
                      cc_or_quant = "cc"),
  eqtl_dir = here::here("raw_data", 
                         "ota_cell_2021",
                        "full_eqtl_results"),
  eqtl_sample_info_file = here::here("raw_data",
                                     "ota_cell_2021",
                                     "E-GEAD-420.sdrf.txt"),
  output_path = here::here("processed_data",
                           "ota_foo_complete",
                           "coloc"),
  foo_genes_to_test = here::here("raw_data",
                                  "GWAS",
                                  "coloc_munged_data",
                                  "foo_ensembl_gene_ids_overlapping_1Mb_window_hit.RDS")
)



# Load pretidied GWAS data ---------------------------------------------------------------

# GWAS


dir.create(args$output_path, recursive = T)

print("Loading GWAS")
gwas_tidy <- fread(args$munged_gwas)



# Find all genes within +/-1Mb of significant GWAS hits -------------------------

print("Extracting genes of interest")

# Use the following code to extract all genes within +/- 1Mb of significant GWAS hits with the following code. In this case we do not require this because we are only looking in specific loci
# ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas_tidy %>%
#   colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
#                                                   CHR_column = "CHR",
#                                                   BP_column = "BP",
#                                                   mart = 38)


ensembl_gene_names_to_test <- readRDS(args$foo_genes_to_test)
# ensembl_gene_names_to_test <- c("LRRK2",   "BST1",    "MCCC1",  "VPS13C",  "DGKQ", "SNCA", "TMEM175") # These have already been processed
#ensembl_gene_names_to_test <- c("MMRN1", "SNCA-AS1") # Additional genes to run


# Load oto et al eqtl snp information ------------------------------------------

# In this dataset:
## Ref allele is the allele in the reference genome - this is coloc Al2
## Alt allele is the effect allele - this is coloc Al1
## There are blood donations from 416 participants. The various blood cells were isolated from these samples so each eqtl has 416 samples.

print("Load eqtl")

## Select eqtl

eqtl_list <- list.files(args$eqtl_dir)
eqtl_list <- eqtl_list[grepl("_nominal.txt", eqtl_list)]

eqtl_sample_no <- 416

# Note: here using the 1000 genomes project east asian population as it gave more results for MAF than other databases (eg gnomad 2. gnomad3 MafH5.gnomAD.v3.1.2.GRCh38 did not have an east asian population)
mafdb1K <- MafDb.1Kgenomes.phase3.GRCh38

# Coloc loop --------------------------------------------------------------


# cl <- parallel::makeCluster(detectCores()) # Detect number of cores by which to parallelise
# doParallel::registerDoParallel(cl) # Set up number of clusters

for(eqtl_number in 1:length(eqtl_list)) {
  
  print(eqtl_number)
  
  eqtl_sample <- eqtl_list[eqtl_number]
  eqtl_name <- eqtl_sample %>% 
    str_replace(., "_nominal.txt", "") %>%  # Some files have been unzipped, so need to replace both options to generate the name and sample number
    str_replace(., "_nominal.txt.gz", "")
  eqtl_path <- str_c(args$eqtl_dir, "/", eqtl_sample)
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
    dplyr::mutate(gene = str_replace(Gene_id, "\\..*$", "")) %>% 
    dplyr::select(-TSS_position, -Number_of_variants_cis)
  
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
    
    # j = 5
    ensembl_gene_id_to_search <-  cis_genes_matched_gwas_eqtl[j]
    
    
    print(str_c(Sys.time(), " - PD_gwas - ", eqtl_name, "_eqtl - ", ensembl_gene_id_to_search))
    
    # Filter eqtls for snps matching a single editing site
    eqtl_tidy_gene_filtered <-  eqtl %>% 
      dplyr::filter(gene == ensembl_gene_id_to_search)    %>% 
      dplyr::filter(nchar(REF) == 1,  # Remove the SNPs in the eqtl that are insertions/deletions (ie not true SNPs)
                    nchar(ALT) == 1) %>% 
      dplyr::mutate(eqtl_dataset = eqtl_name,
                    CHR = str_replace(CHR, "chr", ""), # Convert to ensembl chr format
                    SNP = str_c(CHR, "_", Variant_position_start),
                    N = eqtl_sample_no,
                    beta = `slope(ALT)`) 
    
    if (nrow(eqtl_tidy_gene_filtered) == 0) {
      print(str_c("No QTLs overlapping: ", ensembl_gene_id_to_search))
      next
    }
    
    eqtl_tidy_gene_filtered <- eqtl_tidy_gene_filtered %>%   
      as_granges(seqnames = CHR,    # Make genomic range to get mafs.
                 start = Variant_position_start,
                 end = Variant_position_end) %>% 
      gscores(x = mafdb1K, ranges = ., pop = "EAS_AF") %>% #Use east asian population as the data is from Japanese donors
      as_tibble() %>% 
      drop_na(EAS_AF) %>% 
      dplyr::select(eqtl_dataset,
                    gene,
                    SNP,
                    beta,
                    p.value = nominal_P_value,
                    Al1 = ALT,
                    Al2 = REF,
                    N,
                    maf = EAS_AF,
                    Gene_name) %>% 
      dplyr::filter(!is.na(maf), #Remove NAs and 0 from maf and beta, 
                    maf != 0,
                    beta != 0) %>%
      dplyr::mutate(maf = dplyr::case_when(maf > 0.5 ~ 1-maf, # flip mafs so they are all <0.5
                                            TRUE ~ maf)) %>% 
    colochelpR::check_coloc_data_format(.,
                                        beta_or_pval = "pval",
                                          check_maf = T) %>% 
      dplyr::distinct(SNP, .keep_all = T)
    
    
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



print("Done")
