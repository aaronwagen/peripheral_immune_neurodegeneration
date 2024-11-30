
library(httr)
library(dplyr)
library(tidyverse)



# Choose working location as jane or nemo
location <- "jane"

if (location == "jane") {
    base_dir <- "/Volumes/lab-gandhis/home/users/wagena/periph_immune_neurodegen_publication/"
} else if (location == "nemo") {
    base_dir <- here::here()
    options(bitmapType='cairo') # Add graphics device for nemo
}


NDD_geneID <- readRDS(file = file.path(base_dir, "/derived_data/NDD_ENSG_gene_ids.RDS"))
NDD_geneID_in_aquino <- readRDS(file = file.path(base_dir, "/derived_data/NDD_ENSG_gene_ids_in_aquino.RDS"))
all_gtex_genes <- readRDS(file = file.path(base_dir, "/derived_data/gtex_all_ENSG_ids.RDS"))
all_gtex_genes_noNDD <-   setdiff(all_gtex_genes, NDD_geneID)
brain_gtex_genes <- readRDS(file = file.path(base_dir, "/derived_data/gtex_brain_ENSG_ids.RDS"))
brain_gtex_genes_noNDD <- setdiff(brain_gtex_genes, NDD_geneID)

# Set base URL of GraphQL API endpoint
base_url <- "https://api.genetics.opentargets.org/graphql"
set.seed(613)
n_iterations <- 1000
n_genes_to_check <- length(NDD_geneID_in_aquino)

# GraphQL query for a single gene
l2g_query <- "
query GenePageL2GPipelineQuery($geneId: String!) {
    studiesAndLeadVariantsForGeneByL2G(geneId: $geneId) {
        pval
        yProbaModel
        study {
            studyId
            traitReported
            pubAuthor
            pubDate
            pmid
            nInitial
            nReplication
            hasSumstats
            nCases
        }
        variant {
            rsId
            id
        }
        odds {
            oddsCI
            oddsCILower
            oddsCIUpper
        }
        beta {
            betaCI
            betaCILower
            betaCIUpper
            direction
        }
    }
}"


# Bootstrap all gtex genes ------------------------------------------------

message("Bootstrapping all gtex genes")

bootstrap_chen_aquino_all_gtex_genes <- tibble()

for (iteration_no in 1:n_iterations) {

    message(str_c("Starting iteration: ", iteration_no))
    # List of your gene IDs
    #geneIds <- NDD_geneID_in_aquino# c("ENSG00000145335", " ENSG00000130203", "ENSG00000188906", "ENSG00000165029")
    geneIds <- sample(all_gtex_genes_noNDD, size = n_genes_to_check, replace = F)

    # Add the rest of your gene IDs to the list

    # Initialize an empty list to store results
    all_filtered_data <- list()

    # Loop over each gene ID
    for (geneId in geneIds) {
        # Set variables object of arguments to be passed to endpoint
        variables <- list("geneId" = geneId)

        # Construct POST request body object with query string and variables
        post_body <- list(query = l2g_query, variables = variables)

        # Perform POST request
        r <- POST(url = base_url, body = post_body, encode = 'json')

        # Check for successful response
        if (r$status_code == 200) {
            df <- content(r)
        } else {
            warning("Failed to retrieve data for gene ", geneId, ". Status code: ", r$status_code)
            next  # Skip to the next gene ID
        }

        # Extract the data
        data_list <- df$data$studiesAndLeadVariantsForGeneByL2G

        # Filter and process the data
        filtered_data <- lapply(data_list, function(item) {
            # Check if the item's study pubAuthor matches "Chen MH"
            if (!is.null(item$study$pubAuthor) && item$study$pubAuthor == "Chen MH") {
                item$geneId <- geneId  # Add geneId to the item
                return(item)
            } else {
                return(NULL)  # Return NULL for items that don't match
            }
        })

        # Remove NULL entries from the list
        filtered_data <- Filter(Negate(is.null), filtered_data)

        # Append to the master list
        all_filtered_data <- c(all_filtered_data, filtered_data)
    }

    # Now, convert all_filtered_data to a dataframe
    # iteration_no = 1
    if (length(all_filtered_data) > 0) {
        filtered_df <- do.call(rbind, lapply(all_filtered_data, function(item) {
            # Extract the relevant fields and handle NULLs by replacing them with NA
            data.frame(
                iteration = iteration_no,
                gene = item$geneId,
                pval = ifelse(is.null(item$pval), NA, item$pval),
                studyId = ifelse(is.null(item$study$studyId), NA, item$study$studyId),
                traitReported = ifelse(is.null(item$study$traitReported), NA, item$study$traitReported),
                pubAuthor = ifelse(is.null(item$study$pubAuthor), NA, item$study$pubAuthor),
                nInitial = ifelse(is.null(item$study$nInitial), NA, item$study$nInitial),
                rsId = ifelse(is.null(item$variant$rsId), NA, item$variant$rsId),
                variantId = ifelse(is.null(item$variant$id), NA, item$variant$id),
                betaCI = ifelse(is.null(item$beta$betaCI), NA, item$beta$betaCI),
                direction = ifelse(is.null(item$beta$direction), NA, item$beta$direction),
                stringsAsFactors = FALSE
            )
        }))
    } else {
        filtered_df <- data.frame()  # Empty dataframe if no data matches the filter
    }

    # Display the dataframe
    #print(filtered_df)


    filtered_df_iteration_summary <- filtered_df %>%
        dplyr::filter(traitReported %in% c("White blood cell count",
                                           "Lymphocyte counts",
                                           "Monocyte count",
                                           "Neutrophil count")) %>%
        dplyr::mutate(traitReported = fct_relevel(traitReported, c("Lymphocyte counts",
                                                                   "Monocyte count",
                                                                   "Neutrophil count",
                                                                   "White blood cell count"))) %>%
        dplyr::count(iteration, traitReported) %>% #, studyId) %>%
        pivot_wider(names_from = traitReported,
                    values_from = n)

    bootstrap_chen_aquino_all_gtex_genes  <- bind_rows(bootstrap_chen_aquino_all_gtex_genes,
                                                       filtered_df_iteration_summary)

}



write_delim(bootstrap_chen_aquino_all_gtex_genes ,
            file  = "/Users/wagena/Documents/chen_aquino_bootstrap_all_gtex_genes.csv")




# Bootstrap brain genes ---------------------------------------------------



bootstrap_chen_aquino_brain_gtex_genes <- tibble()

for (iteration_no in 1:n_iterations) {

    message(str_c("Starting iteration: ", iteration_no))

    # List of your gene IDs
    #geneIds <- NDD_geneID_in_aquino# c("ENSG00000145335", " ENSG00000130203", "ENSG00000188906", "ENSG00000165029")
    geneIds <- sample(brain_gtex_genes_noNDD, size = n_genes_to_check, replace = F)

    # Add the rest of your gene IDs to the list

    # Initialize an empty list to store results
    all_filtered_data <- list()

    # Loop over each gene ID
    for (geneId in geneIds) {
        # Set variables object of arguments to be passed to endpoint
        variables <- list("geneId" = geneId)

        # Construct POST request body object with query string and variables
        post_body <- list(query = l2g_query, variables = variables)

        # Perform POST request
        r <- POST(url = base_url, body = post_body, encode = 'json')

        # Check for successful response
        if (r$status_code == 200) {
            df <- content(r)
        } else {
            warning("Failed to retrieve data for gene ", geneId, ". Status code: ", r$status_code)
            next  # Skip to the next gene ID
        }

        # Extract the data
        data_list <- df$data$studiesAndLeadVariantsForGeneByL2G

        # Filter and process the data
        filtered_data <- lapply(data_list, function(item) {
            # Check if the item's study pubAuthor matches "Chen MH"
            if (!is.null(item$study$pubAuthor) && item$study$pubAuthor == "Chen MH") {
                item$geneId <- geneId  # Add geneId to the item
                return(item)
            } else {
                return(NULL)  # Return NULL for items that don't match
            }
        })

        # Remove NULL entries from the list
        filtered_data <- Filter(Negate(is.null), filtered_data)

        # Append to the master list
        all_filtered_data <- c(all_filtered_data, filtered_data)
    }

    # Now, convert all_filtered_data to a dataframe
    # iteration_no = 1
    if (length(all_filtered_data) > 0) {
        filtered_df <- do.call(rbind, lapply(all_filtered_data, function(item) {
            # Extract the relevant fields and handle NULLs by replacing them with NA
            data.frame(
                iteration = iteration_no,
                gene = item$geneId,
                pval = ifelse(is.null(item$pval), NA, item$pval),
                studyId = ifelse(is.null(item$study$studyId), NA, item$study$studyId),
                traitReported = ifelse(is.null(item$study$traitReported), NA, item$study$traitReported),
                pubAuthor = ifelse(is.null(item$study$pubAuthor), NA, item$study$pubAuthor),
                nInitial = ifelse(is.null(item$study$nInitial), NA, item$study$nInitial),
                rsId = ifelse(is.null(item$variant$rsId), NA, item$variant$rsId),
                variantId = ifelse(is.null(item$variant$id), NA, item$variant$id),
                betaCI = ifelse(is.null(item$beta$betaCI), NA, item$beta$betaCI),
                direction = ifelse(is.null(item$beta$direction), NA, item$beta$direction),
                stringsAsFactors = FALSE
            )
        }))
    } else {
        filtered_df <- data.frame()  # Empty dataframe if no data matches the filter
    }

    # Display the dataframe
    #print(filtered_df)


    filtered_df_iteration_summary <- filtered_df %>%
        dplyr::filter(traitReported %in% c("White blood cell count",
                                           "Lymphocyte counts",
                                           "Monocyte count",
                                           "Neutrophil count")) %>%
        dplyr::mutate(traitReported = fct_relevel(traitReported, c("Lymphocyte counts",
                                                                   "Monocyte count",
                                                                   "Neutrophil count",
                                                                   "White blood cell count"))) %>%
        dplyr::count(iteration, traitReported) %>% #, studyId) %>%
        pivot_wider(names_from = traitReported,
                    values_from = n)

    bootstrap_chen_aquino_brain_gtex_genes  <- bind_rows(bootstrap_chen_aquino_brain_gtex_genes,
                                                       filtered_df_iteration_summary)

}



write_delim(bootstrap_chen_aquino_brain_gtex_genes ,
            file  = "/Users/wagena/Documents/chen_aquino_bootstrap_brain_gtex_genes.csv")





