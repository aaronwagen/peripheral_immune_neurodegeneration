# Bootstrap chen asociations
# Aaron Wagen
# A script to bootstrap the number of associations found in Chen compared to the expressed Aquino genes found.

# To run on slurm: 
# sbatch --time=1-00:00:00 --mem=6G --wrap="Rscript r_scripts/opentargets_chen_bootstrap_parallelised.R"

# Load necessary libraries
library(httr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(forcats)
library(doSNOW)
library(progress)
library(logging)

# Choose working location as 'jane' or 'nemo'
location <- "nemo"

if (location == "jane") {
    base_dir <- "/Volumes/lab-gandhis/home/users/wagena/periph_immune_neurodegen_publication/"
} else if (location == "nemo") {
    base_dir <- "/camp/home/wagena/periph_immune_neurodegen_publication/"
    options(bitmapType = 'cairo') # Add graphics device for nemo
}

# Initialize logging
output_dir <- file.path(base_dir, "/derived_data/")
basicConfig()
addHandler(writeToFile, file = file.path(output_dir, "bootstrap_log.txt"))

# Load gene lists
NDD_geneID <- readRDS(file = file.path(base_dir, "/derived_data/NDD_ENSG_gene_ids.RDS"))
NDD_geneID_in_aquino <- readRDS(file = file.path(base_dir, "/derived_data/NDD_ENSG_gene_ids_in_aquino.RDS"))
all_gtex_genes <- readRDS(file = file.path(base_dir, "/derived_data/gtex_all_ENSG_ids.RDS"))
all_gtex_genes_noNDD <- setdiff(all_gtex_genes, NDD_geneID)
brain_gtex_genes <- readRDS(file = file.path(base_dir, "/derived_data/gtex_brain_ENSG_ids.RDS"))
brain_gtex_genes_noNDD <- setdiff(brain_gtex_genes, NDD_geneID)

# Set parameters
set.seed(613)
n_iterations <- 12
n_genes_to_check <- length(NDD_geneID_in_aquino)
traits_of_interest <- c("White blood cell count",
                        "Lymphocyte counts",
                        "Monocyte count",
                        "Neutrophil count")

# Define the function for a single iteration
bootstrap_iteration <- function(iteration_no, gene_pool, n_genes_to_check) {
    # Load necessary packages within the function
    library(httr)
    library(dplyr)
    library(tidyr)
    library(forcats)
    library(logging)
    
    # Define base URL and query
    base_url <- "https://api.genetics.opentargets.org/graphql"
    
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
    
    # Sample genes
    geneIds <- sample(gene_pool, size = n_genes_to_check, replace = FALSE)
    
    # Initialize an empty list to store results
    all_filtered_data <- list()
    
    # Loop over each gene ID
    for (geneId in geneIds) {
        # Set variables object of arguments to be passed to endpoint
        variables <- list("geneId" = geneId)
        
        # Construct POST request body object with query string and variables
        post_body <- list(query = l2g_query, variables = variables)
        
        # Perform POST request
        r <- tryCatch({
            POST(url = base_url, body = post_body, encode = 'json')
        }, error = function(e) {
            logwarn(sprintf("Iteration %d: Failed to retrieve data for gene %s. Error: %s", iteration_no, geneId, e$message))
            return(NULL)
        })
        
        # Check for successful response
        if (is.null(r) || r$status_code != 200) {
            logwarn(sprintf("Iteration %d: Failed to retrieve data for gene %s. Status code: %s", iteration_no, geneId, r$status_code))
            next  # Skip to the next gene ID
        }
        
        df <- content(r)
        
        # Extract the data
        data_list <- df$data$studiesAndLeadVariantsForGeneByL2G
        
        if (is.null(data_list) || length(data_list) == 0) {
            # If data_list is NULL or empty, skip to the next gene
            next
        }
        
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
        # Create an empty dataframe with the same columns
        filtered_df <- data.frame(
            iteration = integer(0),
            gene = character(0),
            pval = numeric(0),
            studyId = character(0),
            traitReported = character(0),
            pubAuthor = character(0),
            nInitial = numeric(0),
            rsId = character(0),
            variantId = character(0),
            betaCI = character(0),
            direction = character(0),
            stringsAsFactors = FALSE
        )
    }
    
    # Process the dataframe
    # Filter only for traits of interest
    filtered_df <- filtered_df %>%
        dplyr::filter(traitReported %in% traits_of_interest)
    
    # Count the occurrences of each trait
    trait_counts <- filtered_df %>%
        dplyr::count(traitReported) %>%
        tidyr::pivot_wider(names_from = traitReported,
                           values_from = n,
                           values_fill = 0)
    
    # Ensure all traits of interest are included in the final output
    for (trait in traits_of_interest) {
        if (!(trait %in% colnames(trait_counts))) {
            trait_counts[[trait]] <- 0
        }
    }
    
    # Add iteration number to the summary
    trait_counts <- trait_counts %>%
        mutate(iteration = iteration_no) %>%
        select(iteration, all_of(traits_of_interest))
    
    # If there was no data, create a row with zeros
    if (nrow(trait_counts) == 0) {
        trait_counts <- data.frame(
            iteration = iteration_no,
            matrix(0, nrow = 1, ncol = length(traits_of_interest)),
            stringsAsFactors = FALSE
        )
        colnames(trait_counts)[-1] <- traits_of_interest
    }
    
    # Return the result
    return(trait_counts)
}

# Function to run bootstrapping
run_bootstrap <- function(gene_pool, output_file) {
    # Number of cores to use
    numCores <- 8  # Adjust this number based on your system and API limitations
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
    
    # Initialize progress bar
    pb <- txtProgressBar(max = n_iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Run the foreach loop with progress bar
    bootstrap_results <- foreach(iteration_no = 1:n_iterations,
                                 .combine = bind_rows,
                                 .options.snow = opts,
                                 .packages = c("httr", "dplyr", "tidyr", "forcats", "logging"),
                                 .export =  c("bootstrap_iteration", "n_genes_to_check", "traits_of_interest")) %dopar% {
        bootstrap_iteration(iteration_no, gene_pool = gene_pool, n_genes_to_check = n_genes_to_check)
    }
    
    # Close progress bar and cluster
    close(pb)
    stopCluster(cl)
    
    # Write the results to a file
    write_csv(bootstrap_results,
                file = file.path(output_dir, output_file))
}

# Run bootstrapping for all GTEx genes
message("Bootstrapping all GTEx genes")
run_bootstrap(all_gtex_genes_noNDD, "chen_aquino_bootstrap_all_gtex_genes_parallel.csv")

# Run bootstrapping for brain GTEx genes
message("Bootstrapping brain GTEx genes")
run_bootstrap(brain_gtex_genes_noNDD, "chen_aquino_bootstrap_brain_gtex_genes_parallel.csv")
