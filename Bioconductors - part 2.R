# Step 1: Install and load the required packages
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# Load the biomaRt package
library(biomaRt)

# Example: Load the ALL dataset if needed
# (Uncomment this if ALL dataset is part of an external package)
# if (!requireNamespace("ALL", quietly = TRUE)) {
#     BiocManager::install("ALL")
# }
# library(ALL)

# Step 2: Connect to Ensembl 75 (hg19 build)
mart <- useMart(host = 'https://feb2014.archive.ensembl.org', 
                biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

# Step 3: Get probeset IDs from the ALL dataset (adjust this step based on your dataset)
# For example, you can extract probeset IDs if the dataset has a featureNames() function:
# all_probeset_ids <- featureNames(ALL)  # Adjust based on your actual dataset

# Here, replace `all_probeset_ids` with the list of actual probeset IDs you want to annotate
all_probeset_ids <- c("202429_at", "202431_s_at", "202433_s_at")  # Example probeset IDs


annotated_probes <- getBM(attributes = c('affy_hg_u133_plus_2', 'ensembl_gene_id'), 
                          filters = 'affy_hg_u133_plus_2',
                          values = all_probeset_ids, 
                          mart = mart)

# Check how many probesets are mapped to more than one gene
probeset_gene_count <- aggregate(ensembl_gene_id ~ affy_hg_u133_plus_2, data = annotated_probes, FUN = length)
probesets_with_multiple_genes <- probeset_gene_count[probeset_gene_count$ensembl_gene_id > 1, ]
num_probesets_with_multiple_genes <- nrow(probesets_with_multiple_genes)

# Output the result
cat("Number of probesets with more than one Ensembl gene ID:", num_probesets_with_multiple_genes, "\n")



## Question 3

mart <- useMart(host = 'https://feb2014.archive.ensembl.org', 
                biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

# Step 3: Get probeset IDs from your dataset
# Example: Use probeset IDs from a dataset like ALL
# all_probeset_ids <- featureNames(ALL)  # Replace with your actual dataset

# Simulating with 1000 example probeset IDs (replace with actual data from your dataset)
all_probeset_ids <- sample(paste0("20", sample(100000:999999, 1000, replace = TRUE), "_at"), 1000)

# Step 4: Annotate probesets with Ensembl gene IDs and chromosome information
# Query for probeset IDs, corresponding Ensembl gene IDs, and chromosome locations
annotated_probes <- getBM(attributes = c('affy_hg_u133_plus_2', 'ensembl_gene_id', 'chromosome_name'), 
                          filters = 'affy_hg_u133_plus_2',
                          values = all_probeset_ids,
                          mart = mart)

# Step 5: Print the first few rows to inspect the data
print(head(annotated_probes))

# Step 6: Filter for genes on autosomes (chromosomes 1 to 22)
# Convert chromosome_name to character to handle any numeric or string inconsistencies
annotated_probes$chromosome_name <- as.character(annotated_probes$chromosome_name)

# Check the unique chromosome names to verify if 1-22 are represented as expected
print(unique(annotated_probes$chromosome_name))

# Filter only for chromosomes 1 to 22 (autosomes)
autosome_genes <- annotated_probes[annotated_probes$chromosome_name %in% as.character(1:22), ]

# Step 7: Count how many unique probesets are annotated with one or more genes on autosomes
unique_probesets_on_autosomes <- length(unique(autosome_genes$affy_hg_u133_plus_2))

# Output the result
cat("Number of probesets annotated with one or more genes on autosomes:", unique_probesets_on_autosomes, "\n")




## Question 4 

# Install and load the minfiData package if not already installed
if (!requireNamespace("minfiData", quietly = TRUE)) {
  BiocManager::install("minfiData")
}

# Load the minfiData package
library(minfiData)

# Load the MsetEx dataset
data(MsetEx)

# Extract the Methylation values
meth_values <- getMeth(MsetEx)

# Extract Methylation values for the sample "5723646052_R04C01"
sample_meth_values <- meth_values[, "5723646052_R04C01"]

# Calculate the mean of the Methylation values across all features for the sample
mean_meth_value <- mean(sample_meth_values, na.rm = TRUE)

# Output the result
cat("Mean Methylation value for sample 5723646052_R04C01:", mean_meth_value, "\n")


## Question 5

# Step 1: Install and load the GEOquery package
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

# Load the GEOquery package
library(GEOquery)

# Step 2: Download the processed data for GSE788
gse <- getGEO("GSE788", GSEMatrix = TRUE)

# Check if multiple versions of the data were downloaded
if (length(gse) > 1) {
  gse <- gse[[1]]  # Select the first set if multiple exist
}

# Step 3: Extract the expression data
expr_data <- exprs(gse)

# Step 4: Extract the expression values for sample GSM9024
sample_gsm9024 <- expr_data[, "GSM9024"]

# Step 5: Calculate the mean expression level for sample GSM9024
mean_expr_gsm9024 <- mean(sample_gsm9024, na.rm = TRUE)

# Output the result
cat("Mean expression level for sample GSM9024:", mean_expr_gsm9024, "\n")


## Question 6 

# Step 1: Install and load the airway package
if (!requireNamespace("airway", quietly = TRUE)) {
  BiocManager::install("airway")
}

# Load the airway package
library(airway)

# Load the airway dataset
data(airway)

# Step 2: Check the structure of the airway dataset
# The dataset is a SummarizedExperiment object, we will access the relevant information
print(airway)

# Extract the read lengths (assuming the dataset contains this information)
# For this example, I'm assuming there is a metadata column with the length of the reads.
# This is just an example to access the metadata (adjust if necessary):
sample_lengths <- rowData(airway)$Length  # Replace 'Length' with the correct metadata column if necessary

# Step 3: Calculate the average length for each sample
average_lengths_per_sample <- colMeans(sample_lengths, na.rm = TRUE)

# Step 4: Calculate the average of the averages across all samples
average_of_average_lengths <- mean(average_lengths_per_sample)

# Output the result
cat("The average of the average lengths across the samples is:", average_of_average_lengths, "\n")

### Question 7

# Load the airway package
library(airway)

# Load the airway dataset
data(airway)

# Step 2: Check the structure of the airway dataset
print(airway)

# Step 3: Access the count data
# The count data is stored in the assay slot of the SummarizedExperiment object
count_data <- assay(airway)

# Step 4: Extract the counts for sample SRR1039512
sample_counts <- count_data[, "SRR1039512"]

# Step 5: Count how many Ensembl genes have a count of 1 or more in SRR1039512
num_genes_with_count_1_or_more <- sum(sample_counts >= 1)

# Output the result
cat("Number of Ensembl genes with a count of 1 or more in sample SRR1039512:", num_genes_with_count_1_or_more, "\n")


## 