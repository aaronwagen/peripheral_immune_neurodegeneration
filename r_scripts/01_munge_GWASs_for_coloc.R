# Run coloc on Immune expression qtls and PD GWAS
# Aaron Wagen
# Jan 2023

# Run this script on a slurm cluster with:
# ml R/4.2.0-foss-2021b
# sbatch --time=3-00:00:00 --output logs/coloc_gwas_eqtl.log --open-mode=append --mem=128G --wrap="Rscript r_scripts/run_coloc_gwas_immune.R"


# This script prepares the following GWAS's for use with coloc:
# Nalls 2019 ex 23&me data
# Sakaue and Kanai 2021, Nature Genetics
# Foo et al, 2020, JAMA neurology
# Note the Ota eqtl:
# - do not have ChrX so this will be excluded from the GWASs
# - have chromosomes with 'chr' included, though this is removed in the 02_run_coloc_gwas_immune.R code. I will add the chr for the purpose of the liftover from GRch37 to GRCh38, and then remove again to run coloc.
# - the SNP I have used in the run coloc is location based in the format chr_bp, ie: 2_2039422

# GWAS datasets:

## 1) Sakaue and Kanai - Nature Genetics 2021, A cross population atlas of genetic associations for 220 human phenotypes. 
## https://www.nature.com/articles/s41588-021-00931-x#Sec27
## This includes 340 PD cases, and 175788 controls.
## Data downloaded from: https://humandbs.biosciencedbc.jp/en/hum0197-v3-220
## Dictionary file at: https://humandbs.biosciencedbc.jp/files/hum0197.org/Dictfile_BBJ.html


## 2) Foo, Chew, Chung et al - JAMA Neurology 2020, Identification of Risk Loci for Parkinson Disease in Asians and Comparison of Risk Between Asians and Europeans
## https://jamanetwork.com/journals/jamaneurology/fullarticle/2764340
## Columns: CHR, BP, SNP (rs number:location:A1:A2), A1, A2, BETA, P SE, MarkerName (rs number)
## This includes 6724 PD cases, and 24581 controls, 5,843,213 snps
## File supplied by authors via C Blauwendraat.


## 3) Rizig Bandres-Ciga et al, medRxiv, 2023, Genome-wide Association Identifies Novel Etiological Insights Associated with Parkinson’s Disease in African and African Admixed Populations
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10312852/
## Columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, variant_id, ref_allele, direction, HetISq, HetChiSq, HetDf, HetPVal, rsid
## This includes 1488 PD cases, and 196,430 controls
## File supplied by authors via C Blauwendraat.

# Load libraries -------------------------------------------------------------------

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
library(MafDb.1Kgenomes.phase3.GRCh38) # If using GRCh38, if not then use GRCh37.
library(rutils)
library(GenomicScores)

# Note: here using the 1000 genomes project east asian population as it gave more results for MAF than other databases (eg gnomad 2. gnomad3 MafH5.gnomAD.v3.1.2.GRCh38 did not have an east asian population)
mafdb1K <- MafDb.1Kgenomes.phase3.GRCh38




# Define paths ------------------------------------------------------------


# Choose working location as jane or nemo
# location <- "nemo" 
# location <- "jane"

if (location == "jane") {
  base_dir <- "/Volumes/lab-gandhis/home/users/wagena/SNCA"
} else if (location == "nemo") {
  base_dir <- here::here()
}



here::i_am("r_scripts/01_munge_GWASs_for_coloc.R")

args <- list(
  nalls2019ex23andme_gwas = c(path = "/camp/lab/gandhis/home/shared/LDSC/GWAS/PD2019_meta5_ex23andMe/PD2019_ex23andMe",
                 n_cases = 33674,
                 n_controls = 449056,
                 n_total = 482730,
                 prop_cases = 0.074989,
                 cc_or_quant = "cc"),
  sakaue2021_gwas = c(path = file.path(base_dir,
                                       "raw_data",
                                       "GWAS",
                                       "Sakaue_Kanai_Japanese_brain_bank_GWAS_220_traits_PD",
                                       "GWASsummary_Parkinsons_Disease_Japanese_SakaueKanai2020.auto.txt.gz"),
                  n_cases = 340,
                  n_controls = 175788,
                  n_total = 176128,
                  prop_cases = 0.001930,
                  cc_or_quant = "cc"),
  foo2020_gwas = c(path = file.path(base_dir,
                                    "raw_data",
                                    "GWAS",
                                    "foo_chew_etal_JAMA_east_asian_PD_GWAS_sumstats.txt"),
                   n_cases = 6724,
                   n_controls = 24581,
                   n_total = 31305,
                   prop_cases = 0.214790,
                   cc_or_quant = "cc"),
  rizig2023_gwas = c(path = file.path(base_dir,
                                      "raw_data",
                                      "GWAS",
                                      "Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38.txt"),
                     n_cases = 1488,
                     n_controls = 196430,
                     n_total = 197918,
                     prop_cases = 0.007518265,
                     cc_or_quant = "cc"),
  jansen2021_gwas = c(path = file.path(base_dir,
                                       "raw_data",
                                       "GWAS",
                                       "AD_sumstats_Jansenetal_2019sept.txt.gz"),
                      n_cases = 71880, # Note that with this dataset numbers are given for each SNP, so don't need to use these numbers.
                      n_controls = 383378,
                      n_total = 455258,
                      prop_cases = 0.157888494,
                      cc_or_quant = "cc"),
  shigemizu2021_gwas = c(path = file.path(base_dir,
                                          "raw_data",
                                          "GWAS",
                                          "shigemizu_2021_japanese_AD_GWAS.txt"),
                         n_cases = 3962,
                         n_controls = 4074,
                         n_total = 8037,
                         prop_cases = 0.4930313589,
                         cc_or_quant = "cc"),
  hg19tohg38_liftover = "/nemo/lab/gandhis/home/users/wagena/references/GRCh38/tools/liftover/hg19ToHg38.over.chain",
  example_ota_eqtl = file.path(base_dir, "raw_data/ota_cell_2021/CD16p_Mono_nominal.txt"),
  gwas_results_path = file.path(base_dir, "raw_data/GWAS/coloc_munged_data")
  )


# Create GWAS summary dataframe from list above
gwas_summary <- as.data.frame(t(do.call(cbind, args))) %>% 
  rownames_to_column(var = "gwas") %>% 
  dplyr::mutate(gwas = str_replace(gwas, "_gwas", ""))


# Define functions --------------------------------------------------------


rename_gwas_columns <- function(gwas_df) {
  #' Function that renames GWAS columns to all be in the same format, allowing other functions to run. Note, it will not check whether Al1 and Al2 are the correct reference
  #' @param gwas_df dataframe from a GWAS with columns to rename. Specifically will use the following labels: beta, p.value, Al1, Al2, 
  #' @
  #' @return dataframe with renamed columns

  
  gwas_df <- gwas_df %>% 
    dplyr::rename(beta = any_of(c("BETA", "beta", "b")),
                  p.value = any_of(c( "p.value", "p", "P", "P.VALUE", "p_value")),
                  se = any_of(c( "SE", "se", "standard.error", "standard_error")),
                  Al1 = any_of(c( "A1", "Al1", "Allele1")),
                  Al2 = any_of(c( "A2", "Al2", "Allele2"))
    )
  
  return(gwas_df)
}



generate_location_based_SNP_column <- function(gwas_df) {
  #' Function that takes a dataframe with columns for CHR and BP, and generates a location based SNP column based on the format CHR_BP
  #' @param gwas_df dataframe from a GWAS containing the columns CHR and BP
  #' @return dataframe with the combined column in the 'SNP' column

  gwas_df <- gwas_df %>% 
    dplyr::mutate(SNP = str_c(CHR, "_", BP))
  
  return(gwas_df)
                  
}


useful_gwas_numbers <- function(gwas_df, maf_present = T) {
  #' Function that reports useful qc numbers while munging a GWAS
  #' @param gwas_df dataframe from a GWAS containing the columns CHR and BP
  #' @param maf_present, TRUE or FALSE value whether a maf column is present in the dataframe
  #' @return List of useful qc numbers
  

    
    print(stringr::str_c("Total number of SNPs: ",
                         nrow(gwas_df)))
    
    
    print(stringr::str_c("Number of multiallelic SNPs: ",
                         duplicated(gwas_df$SNP) %>%
                           sum()))
    
    print(stringr::str_c("Number of SNPs with beta == 0: ",
                         gwas_df %>%
                           dplyr::filter(beta == 0) %>%
                           nrow()))
    
    print(stringr::str_c("Number of SNPs with multiple nucleotides (Indels): ",
                         sum(nchar(gwas_df$Al1)>1, nchar(gwas_df$Al2)>1)))
                         
    
    if (maf_present == TRUE) {
    
    print(stringr::str_c("Number of SNPs with MAF = NA: ",
                         gwas_df %>%
                           dplyr::filter(is.na(maf)) %>%
                           nrow()))
    
    print(stringr::str_c("Number of SNPs with MAF > 0.5: ",
                         gwas_df %>%
                           dplyr::filter(maf > 0.5) %>%
                           nrow()))
  }
}
  



remove_snps_beta_0 <- function(gwas_df) {
  #' Function that removes values where beta = 0 or NA
  #' @param gwas_df dataframe from a GWAS containing the columns beta
  #' @return dataframe with SNPs removed where beta = 0

  gwas_df <- gwas_df %>%
    dplyr::filter(beta != 0) 
  
  return(gwas_df)
}



remove_duplicated_snps <- function(gwas_df) {
  #' Function that removes duplicated snps by chosing the one with the lowest p value
  #' @param gwas_df dataframe from a GWAS, with the column SNP and p value
  #' @return dataframe with duplicated SNPs removed

# Find duplicated SNPs
duplicated_snps <- gwas_df %>%
  dplyr::filter(duplicated(SNP)) %>%
  dplyr::select(SNP) %>%
  unlist()

# Choose duplicate with lowest p value
filtered_duplicates <- gwas_df %>%
  dplyr::filter(SNP %in% duplicated_snps) %>%
  dplyr::group_by(SNP) %>%
  dplyr::slice_min(order_by = p.value, n=1, with_ties = F) %>%
  dplyr::ungroup()

# Add in filtered duplicates to non-duplicated snps
gwas_tidy <- gwas_df %>%
  dplyr::filter(!SNP %in% duplicated_snps) %>% # Remove all duplicates previously identified
  dplyr::bind_rows(filtered_duplicates)  # Add the rows of filtered duplicates


return(gwas_tidy)

}


change_indels <- function(gwas_df) {
  #' Function that recodes indels (SNPs with multiple nucleotides in Al1 or Al2) as Insertions and Deletions
  #' @param gwas_df dataframe from a GWAS, with the column Al1 and Al2 for Allele 1 and 2
  #' @return dataframe with duplicated SNPs removed
  
  gwas_df$Al2[nchar(gwas_df$Al1)>1] = "D"
  gwas_df$Al1[nchar(gwas_df$Al1)>1] = "I"
  gwas_df$Al1[nchar(gwas_df$Al2)>1] = "D"
  gwas_df$Al2[nchar(gwas_df$Al2)>1] = "I"
  
  return(data)
  
}

retrieve_mafs <- function(gwas_df, population = c("AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF")) {
  #' Function that retrieves the minor allele frequence for a SNP based on the 1000 genomes project phase 3.
  #' @param gwas_df dataframe from a GWAS, with the columns CHR and BP from the genome build matching the one used in the maf database loaded
  #' @param population which of the [geographic populations](https://www.ensembl.org/Help/Faq?id=532) including: African (AFR_AF);  Ad Mixed American (AMR_AF); East Asian (EAS_AF); European (EUR_AF); and south asian (SAS_AF)
  #' @return dataframe with mafs, and those rows that do not have a maf removed.
  gwas <- gwas %>% 
  as_granges(seqnames = CHR,    # Make genomic range to get mafs.
             start = BP,
             end = BP) %>% 
  gscores(x = mafdb1K, ranges = ., pop = population) %>% #Use east asian population as the data is from Japanese donors
  as_tibble() %>% 
  drop_na(EAS_AF) %>% 
  dplyr::rename(maf = EAS_AF,
                CHR = seqnames,
                BP = start) %>% 
  dplyr::select(!c(end, width, strand))
}
  
  
  

# Nalls 2019 GWAS  ---------------------------------------------------------------


# Nalls GWAS: Allele 1 in the GWAS is the effect allele, and allele 2 is the reference (GRCh38) allele - matching A1 and A2 respectively
# This is mapped to GRCh38 reference.

# Note GWAS has been tidied by: 
#   Removing SNPs with beta 0; 
#   remove duplicate SNP entries by choosing the SNP with the lowest p value; 
#   Switching all mafs removing all maf's > 0.5 by subtracting these from 1; 
#   Adding a varbeta (variation in beta) variable for each SNP - a marker of precision of the estimate - using the colochelpR package
#   Adding in number variables, including total number of participants in the GWAS and proportion of cases verse controls.



# Code to tidy GWAS
gwas_name = "nalls2019ex23andme"
gwas_path <- gwas_summary[gwas_summary$gwas == gwas_name, "path"]

gwas <- fread(gwas_path)
gwas_orig <- gwas
gwas <- gwas_orig %>%
  dplyr::mutate(GWAS = "PD2019_ex23andMe",
                N = args$PD_gwas_details$n_total,
                N_propor = args$PD_gwas_details$prop_cases,
                gwas_locus = str_c(CHR,"_", BP),
                CHR = as.factor(CHR)) %>%
  dplyr::select(GWAS,
                CHR,
                BP,
                gwas_locus,
                SNP,
                beta = b,
                se,
                p.value = p,
                Al1 = A1,
                Al2 = A2,
                maf = freq,
                N,
                N_propor)

useful_gwas_numbers(gwas)

# Remove NAs and 0 from maf and beta, and flip mafs to they are all <0.5
gwas_tidy <-
  gwas %>%
  dplyr::filter(!is.na(maf),
                beta != 0) %>%
  dplyr::mutate(
    maf =
      dplyr::case_when(
        maf > 0.5 ~ 1-maf,
        TRUE ~ maf
      )
  )


# Find duplicated SNPs
duplicated_snps <- gwas %>%
  dplyr::filter(duplicated(SNP)) %>%
  dplyr::select(SNP) %>%
  unlist

# Choose duplicate with lowest p value
filtered_duplicates <- gwas_tidy %>%
  dplyr::filter(SNP %in% duplicated_snps) %>%
  dplyr::group_by(gwas_locus) %>%
  dplyr::slice_min(order_by = p.value, n=1, with_ties = T) %>%
  dplyr::ungroup()


gwas_tidy <- gwas_tidy %>%
  dplyr::filter(!SNP %in% duplicated_snps) %>% # Remove all duplicates previously identified
  dplyr::bind_rows(filtered_duplicates) %>%  # Add the rows of filtered duplicates
  dplyr::arrange(CHR, BP) %>%
  colochelpR::get_varbeta(.) %>%
  colochelpR::check_coloc_data_format(beta_or_pval = "beta", check_maf = F) # Check columns are the correct format


useful_gwas_numbers(gwas_tidy)

write_delim(gwas_tidy,
            file = str_c(args$gwas_results_path,
                         "/PD2019_ex23andMe_tidy_varbeta.txt"))



# Find all genes within +/-1Mb of significant GWAS hits
# I have pre-saved the output as the ensembl portal often fails.

ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas_tidy %>%
  dplyr::mutate(CHR = as.integer(CHR)) %>%
  colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                                  CHR_column = "CHR",
                                                  BP_column = "BP",
                                                  mart = 38)

saveRDS(ensembl_gene_ids_overlapping_1Mb_window_hit,
        file = args$ensembl_gene_ids_overlapping_1Mb_window_hit_path)  


# Foo 2020 GWAS -----------------------------------------------------------




## 2) Foo, Chew, Chung et al - JAMA Neurology 2020, Identification of Risk Loci for Parkinson Disease in Asians and Comparison of Risk Between Asians and Europeans
## https://jamanetwork.com/journals/jamaneurology/fullarticle/2764340
## Columns: CHR, BP, SNP (rs number:location:A1:A2), A1, A2, BETA, P SE, MarkerName (rs number)
## This includes 6724 PD cases, and 24581 controls, 5,843,213 snps
## File supplied by authors via C Blauwendraat.
## Foo GWAS: Allele 2 is the reference (GRCh37) allele, Allele 1 is the effect allele, respectively matching A1 and A2 required for coloc.
## This is mapped to GRCh37 reference (note the saved tidy_varbeta file has been mapped to GRCh38)

# Note GWAS has been tidied by: 
#   Lifting over the GRCh37 locations to GRCh38
#   Removing SNPs with beta 0; 
#   Switching all mafs removing all maf's > 0.5 by subtracting these from 1; 
#   Adding a varbeta (variation in beta) variable for each SNP - a marker of precision of the estimate - using the colochelpR package
#   Adding in number variables, including total number of participants in the GWAS and proportion of cases verse controls.




gwas_name = "foo2020"
gwas_path <- gwas_summary[gwas_summary$gwas == gwas_name, "path"]

foo2020 <- fread(gwas_path)

gwas <- foo2020 %>% 
  dplyr::mutate(GWAS = gwas_name,
                CHR = str_c("chr", CHR)) %>%  # Add in "chr" string to use for liftover function, and to match the eqtl)
  dplyr::select(-SNP,
                rs_number = MarkerName) %>% 
  rename_gwas_columns(.)

              
  
# If GRCH37 use liftover command
gwas <- rutils::liftover_coord(df = gwas,  # Liftover from GRCh37 to 38, note this function removes the "chr" prefix, which is what we want to match our code for coloc
                               path_to_chain = args$hg19tohg38_liftover)

# If using location based SNP column for coloc, generate it using the command below from CHR and BP columns
gwas <- generate_location_based_SNP_column(gwas)

# Generate varbeta
gwas <- colochelpR::get_varbeta(gwas) %>% 
  drop_na(varbeta) # Remove SNPs that do not have a SE/varbeta

gwas <- retrieve_mafs(gwas, population = "EAS_AF")


useful_gwas_numbers(gwas, maf_present = TRUE)

# Remove NAs and 0 from maf and beta
gwas <- remove_snps_beta_0(gwas) 

# Remove multiallelic snps by taking the snp with the lowest p value
gwas <- remove_duplicated_snps(gwas)

colochelpR::check_coloc_data_format(gwas, beta_or_pval = "beta", check_maf = T) # Check columns are the correct format


write_delim(gwas,
            file = str_c(args$gwas_results_path,
                         "/foo2020_gwas_tidy_varbeta.txt"))




# Find all genes within +/-1Mb of significant GWAS hits
# I have pre-saved the output as the ensembl portal often fails.

foo_ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas %>%
  dplyr::mutate(CHR = as.integer(CHR)) %>%
  colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                                  CHR_column = "CHR",
                                                  BP_column = "BP",
                                                  mart = 38)

saveRDS(foo_ensembl_gene_ids_overlapping_1Mb_window_hit,
        file = str_c(args$gwas_results_path,
                     "/foo_ensembl_gene_ids_overlapping_1Mb_window_hit.RDS"))

# Rizig 2023 African ------------------------------------------------------

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10312852/
## Columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, variant_id, ref_allele, direction, HetISq, HetChiSq, HetDf, HetPVal, rsid

## Columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, variant_id, ref_allele, direction, HetISq, HetChiSq, HetDf, HetPVal, rsid
## This includes 1488 PD cases, and 196,430 controls
## File supplied by authors via C Blauwendraat.
## There are 3 columns relating to alleles: effect_allele, other_allele, and ref_allele. Coloc requires Allele 2 as the reference allele (in GRCh38), and Allele 1 as the effect allele. 
##    This means that: where the ref_allele == other_allele, ref_allele=Al2, effect_allele=Al1
##                    where the ref_allele == effect_allele, ref_allele=Al2, other_allele=Al1, beta = -beta, (and effect_allele_frequency = 1-effect_allele_frequency)
## This is mapped to GRCh38 reference

# Note GWAS has been tidied by: 
#   Columns renamed
#   Assign effect and ref alleles as above
#   Removing SNPs with beta 0; 
#   Switching all mafs that are > 0.5 to below 0.5 by subtracting these from 1; 
#   Adding a varbeta (variation in beta) variable for each SNP - a marker of precision of the estimate - using the colochelpR package
#   Adding in number variables, including total number of participants in the GWAS and proportion of cases verse controls.



gwas_name = "rizig2023"
gwas_path <- gwas_summary[gwas_summary$gwas == gwas_name, "path"]

rizig2023 <- fread(gwas_path)


gwas <- rizig2023 %>% 
  mutate(GWAS = gwas_name,
         Al1 = case_when(
           ref_allele == other_allele ~ effect_allele,
           ref_allele == effect_allele ~ other_allele,
           TRUE ~ NA_character_  # Handle the case where neither condition is true
         ),
         Al2 = case_when(
           ref_allele == other_allele ~ ref_allele,
           ref_allele == effect_allele ~ ref_allele,
           TRUE ~ NA_character_  # Handle the case where neither condition is true
         ),
         beta = case_when(
           ref_allele == effect_allele ~ -beta,
           TRUE ~ beta  # Keep original beta otherwise
         ),
         effect_allele_frequency = case_when(
           ref_allele == effect_allele ~ 1 - effect_allele_frequency,
           TRUE ~ effect_allele_frequency  # Keep original frequency otherwise
         )) %>% 
  rename_gwas_columns(.) %>% 
  dplyr::select(GWAS,
                CHR = chromosome,
                BP = base_pair_location,
                SNP = rsid, # Use rs numbers for this to ensure the ref/alt alleles are correct
                beta,
                se,
                p.value,
                Al1,
                Al2,
                maf = effect_allele_frequency)  %>% 
  dplyr::mutate(maf = case_when(
    maf <=0.5 ~ maf,
    maf >0.5 ~ 1- maf
  )) %>% 
  colochelpR::get_varbeta(.) %>% 
  drop_na(varbeta) %>% 
  remove_snps_beta_0(.) %>%  
  remove_duplicated_snps(.)


useful_gwas_numbers(gwas, maf_present = TRUE)
colochelpR::check_coloc_data_format(gwas, beta_or_pval = "beta", check_maf = T) # Check columns are the correct format



write_delim(gwas,
            file = str_c(args$gwas_results_path,
                         "/rizig2023_gwas_tidy_varbeta.txt"))


# Find all genes within +/-1Mb of significant GWAS hits
# I have pre-saved the output as the ensembl portal often fails.

rizig_ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas %>%
  dplyr::mutate(CHR = as.integer(CHR)) %>%
  colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                                  CHR_column = "CHR",
                                                  BP_column = "BP",
                                                  mart = 38)

saveRDS(rizig_ensembl_gene_ids_overlapping_1Mb_window_hit,
        file = str_c(args$gwas_results_path,
                     "/rizig_ensembl_gene_ids_overlapping_1Mb_window_hit.RDS"))





# Jansen 2019 AD ----------------------------------------------------------

# https://www.nature.com/articles/s41588-018-0311-9
# AD GWAS downloaded from: https://cncr.nl/research/summary_statistics/
# Note: this is mapped using GRCh37. However, because rs numbers and mafs are listed, we do not need to do a liftover to GRCh38.

# As per the readme on the website the columns are as follows:

# uniqID.a1a2 = uniq ID per variant, format -> CHR:BP_A1_A2
# CHR = chromosome
# BP = base pair position
# A1 = allele 1 (effect allele)
# A2 = allele 2 (non-effect allele) - the gene in GRCh37 reference genome
# SNP = rsID (if available)
# Z = Z-statistic
# P = p-value
# Nsum = sample size by simply summing n per cohort
# Neff = effecting sample size
# dir = directions of effect per cohort (order: ADSP, IGAP, UKB, PGC-ALZ)
# EAF = allele frequency of allele 1
# BETA = effect size
# SE = standard error

# Note that with this dataset numbers are given for each SNP, so don't need to use the total numbers included in the study. 
# I will use Neff as the effective N, which takes into account potential similarities between study participants moreso than Nsum.
# I will derive the proportion of SNPs affected from the overall numbers in the study and assume this is stable across the SNPs.


gwas_name = "jansen2021"
gwas_path <- gwas_summary[gwas_summary$gwas == gwas_name, "path"]

jansen2021 <- fread(gwas_path)

gwas_tidy <- jansen2021 %>% 
  as_tibble() %>% 
  rename_gwas_columns(.) %>% 
  dplyr::rename(maf = EAF) %>% 
  dplyr::mutate(GWAS = gwas_name,
                prop_cases = args$jansen2021_gwas[["prop_cases"]]) %>% 
  dplyr::select(GWAS,
                SNP,
                beta,
                se,
                p.value,
                Al1,  
                Al2,
                maf,
                N = Neff,
                prop_cases,
                CHR,
                BP_grch37 = BP) %>% 
  remove_snps_beta_0(.) %>% 
  remove_duplicated_snps(.) %>% 
  dplyr::mutate(maf = case_when(
    maf <=0.5 ~ maf,
    maf >0.5 ~ 1- maf)) %>% 
  get_varbeta() %>% 
  dplyr::filter(!is.na(maf),
                !is.na(varbeta))




useful_gwas_numbers(gwas_tidy, maf_present = TRUE)
colochelpR::check_coloc_data_format(gwas_tidy, beta_or_pval = "beta", check_maf = T) # Check columns are the correct format


write_delim(gwas_tidy,
            file = str_c(args$gwas_results_path,
                         "/jansen2019_gwas_tidy_varbeta.txt"))




# Find all genes within +/-1Mb of significant GWAS hits
# I have pre-saved the output as the ensembl portal often fails.

jansen_ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas_tidy %>%
  dplyr::mutate(CHR = as.integer(CHR)) %>%
  colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                                  CHR_column = "CHR",
                                                  BP_column = "BP_grch37",
                                                  mart = 37)

saveRDS(jansen_ensembl_gene_ids_overlapping_1Mb_window_hit,
        file = str_c(args$gwas_results_path,
                     "/jansen_ensembl_gene_ids_overlapping_1Mb_window_hit.RDS"))




# Shigemizu 2021 Asian AD GWAS -------------------------------------------------


# Shigemizu_2021_japanese_AD_GWAS.txt
# - From paper `Ethnic and trans-ethnic genome-wide association studies identify new loci influencing Japanese Alzheimer’s disease risk`
# - https://www.nature.com/articles/s41398-021-01272-3#Sec19
# - The summary stats were shared with me by the author on 18th May 2024
# - This includes 3962 cases and 4074 controls
# - It is aligned to GRCh37
# The columns in the raw data are: 
#   - CHR; SNP; BP (referencing GRCh37) 
# - A1 (reference allele in GRCh37); this will become the coloc AL2
# - A2 effect allele; this will become the new coloc Al1
# - NMISS - number of cases and controls
# - NMISS_A - number of affected; NMISS_U - number of unaffected/ controls
# - MAF_A - maf of effect allele/cases;  MAF_U - maf of control allele/controls
# - OR - odds ratio; SE - standard error of odds ratio
# - L95 - lower limit of 95% confidence interval; U95 -upper limit of 95% confidence interval
# - STAT
# - P - p value
# - Info NCGG; Info Niigata

# For coloc we need the following columns:
# - gwas_name
# - SNP - will use rs number
# - We don't have a beta value so will use maf and p value as metrics for coloc. We will use the raw p value (not a corrected p value)
# - Al1 - the effect allele, currently A2
# - Al2 - the alternate allele, currently A1
# - maf - this is the maf of the affect allele, MAF_A
# - The N and proportion of cases will come directly from the sumstats: N = NMISS, prop_cases = NMISS_A/NMISS

gwas_name = "shigemizu2021"
gwas_path <- gwas_summary[gwas_summary$gwas == gwas_name, "path"]

shigemizu2021 <- fread(gwas_path)

gwas_tidy <- shigemizu2021 %>% 
  dplyr::rename(Al1 = A2,
                Al2 = A1) %>% 
  rename_gwas_columns(.) %>% 
  dplyr::mutate(GWAS = gwas_name,
                prop_cases = NMISS_A/NMISS) %>% 
  dplyr::select(GWAS,
                SNP,
                p.value,
                Al1,  
                Al2,
                maf = MAF_A,
                N = NMISS,
                prop_cases,
                CHR,
                BP_GRCh37 = BP) %>% 
  dplyr::mutate(maf = case_when(
    maf <=0.5 ~ maf,
    maf >0.5 ~ 1- maf)) 
  
colochelpR::check_coloc_data_format(gwas_tidy, beta_or_pval = "pval", check_maf = T)

write_delim(gwas_tidy,
            file = str_c(args$gwas_results_path,
                         "/shigemizu2021_gwas_tidy_maf.txt"))



# Find all genes within +/-1Mb of significant GWAS hits
# I have pre-saved the output as the ensembl portal often fails.

shigemuzu_ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas_tidy %>%
  dplyr::mutate(CHR = as.integer(CHR)) %>%
  colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
                                                  CHR_column = "CHR",
                                                  BP_column = "BP_GRCh37",
                                                  mart = 37)

saveRDS(shigemuzu_ensembl_gene_ids_overlapping_1Mb_window_hit,
        file = str_c(args$gwas_results_path,
                     "/shigemuzu_ensembl_gene_ids_overlapping_1Mb_window_hit.RDS"))
