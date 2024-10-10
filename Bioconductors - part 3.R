##Question 1 

# Step 1: Install and load necessary packages
if (!requireNamespace("yeastRNASeq", quietly = TRUE)) {
  BiocManager::install("yeastRNASeq")
}
if (!requireNamespace("ShortRead", quietly = TRUE)) {
  BiocManager::install("ShortRead")
}

# Load the required libraries
library(yeastRNASeq)
library(ShortRead)

# Step 2: Get the FASTQ file path
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")

# Step 3: Read the FASTQ file
fastq_data <- readFastq(fastqFilePath)

# Step 4: Extract the sequences from the FASTQ file
sequences <- sread(fastq_data)

# Step 5: Extract the 5th base from each sequence
fifth_base <- subseq(sequences, start = 5, end = 5)

# Step 6: Calculate the fraction of reads with an "A" at the 5th position
num_reads_with_A_in_5th_base <- sum(fifth_base == "A")
total_reads <- length(fifth_base)
fraction_with_A <- num_reads_with_A_in_5th_base / total_reads

# Output the result
cat("Fraction of reads with an A in the 5th base:", fraction_with_A, "\n")


## Question 2 assessing quality 

quality_scores <- quality(fastq_data)  # This retrieves quality scores

# Step 5: Extract the 5th base quality scores from each read
fifth_base_quality <- as(quality_scores, "matrix")[, 5]

# Step 6: Calculate the average numeric quality score for the 5th base
average_quality_5th_base <- mean(fifth_base_quality, na.rm = TRUE)

# Output the result
cat("Average quality score for the 5th base:", average_quality_5th_base, "\n")


## Question 3 

# Step 1: Install and load necessary packages
if (!requireNamespace("Rsamtools", quietly = TRUE)) {
  BiocManager::install("Rsamtools")
}

if (!requireNamespace("leeBamViews", quietly = TRUE)) {
  BiocManager::install("leeBamViews")
}

# Load the required libraries
library(leeBamViews)
library(Rsamtools)

# Step 2: Get the BAM file path
bamFilePath <- system.file("bam", "isowt5_13e.bam", package = "leeBamViews")

# Step 3: Define the region of interest (chromosome 13, positions 800,000 to 801,000)
region <- GRanges("Scchr13", IRanges(start = 800000, end = 801000))

param <- ScanBamParam(which = region, what = c("pos", "seq"))

# Step 5: Scan the BAM file for reads in the specified region
bam_data <- scanBam(bamFilePath, param = param)


# Step 5: Extract the read positions
read_positions <- bam_data[[1]]$pos  # Extract the positions of reads

# Step 6: Identify duplicated read positions
duplicated_positions <- read_positions[duplicated(read_positions)]

# Step 7: Count how many reads are duplicated by position
num_duplicated_reads <- length(duplicated_positions)

# Output the result
cat("Number of reads duplicated by position in the interval 800,000 to 801,000 on chromosome 13:", num_duplicated_reads, "\n")


## Question 4 


# Load the required libraries
library(leeBamViews)
library(Rsamtools)

# Step 2: Get all BAM file paths
bpaths <- list.files(system.file("bam", package = "leeBamViews"), pattern = "bam$", full.names = TRUE)

# Step 3: Define the region of interest (Scchr13:807762-808068)
region_of_interest <- GRanges("Scchr13", IRanges(start = 807762, end = 808068))

# Step 4: Define the BAM scanning parameters to extract 'pos' (positions of reads)
param <- ScanBamParam(which = region_of_interest, what = c("pos"))

# Step 5: Initialize a variable to hold the total count of reads
total_reads <- 0

# Step 6: Loop over each BAM file and count the reads in the region
for (bamFile in bpaths) {
  # Scan the BAM file for reads in the specified region
  bam_data <- scanBam(bamFile, param = param)
  length(bam_data)
  
  # Extract the read positions
  read_positions <- bam_data[[1]]$pos
 # print(read_positions)
  # Count the number of reads in this sample
  num_reads_in_sample <- length(read_positions)
  print("reads in sample")
  print(num_reads_in_sample)
  
  # Add to the total number of reads
  total_reads <- total_reads + num_reads_in_sample
  #print(total_reads)
}

# Step 7: Calculate the average number of reads across the 8 samples
average_reads <- total_reads / length(bpaths)

# Output the result
cat("Average number of reads in the region Scchr13:807762-808068 across the 8 samples:", average_reads, "\n")


## Question 5

#  Install and load the necessary packages
if (!requireNamespace("oligo", quietly = TRUE)) {
  BiocManager::install("oligo")
}

# Load  package
library(oligo)
library(GEOquery)

getGEOSuppFiles("GSE38792")
list.files("GSE38792")

untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL")

celfiles <- list.files("GSE38792/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)

eset <- rma(rawData)
# Extract the group variable to identify control samples
# Assuming 'group' is a column in the pData (phenotypic data) of the ExpressionSet
#group_info <- pData(eset)$group

# Step 1: Define the group variable manually based on the sample names
group_info <- ifelse(grepl("Control", rownames(pData(eset))), "control", "OSA")

# Step 2: Extract the indices of the control group samples
control_indices <- which(group_info == "control")

# Step 3: Extract the expression values for the probeset "8149273"
probeset_expression <- exprs(eset)["8149273", ]

# Step 4: Calculate the average expression for the control group
control_expression <- probeset_expression[control_indices]
average_control_expression <- mean(control_expression, na.rm = TRUE)

# Output the result
cat("Average expression for the probeset '8149273' in the control group:", average_control_expression, "\n")


## Question 6

# Step 1: Install and load the necessary packages
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

# Load the limma package
library(limma)

# Step 2: Define the group variable (based on previous step)
group_info <- ifelse(grepl("Control", rownames(pData(eset))), "control", "OSA")

# Step 3: Create the design matrix
design <- model.matrix(~ 0 + group_info)
colnames(design) <- c("Control", "OSA")

# Step 4: Fit the linear model using lmFit
fit <- lmFit(exprs(eset), design)

# Step 5: Create a contrast matrix for OSA vs Control
contrast_matrix <- makeContrasts(OSA_vs_Control = OSA - Control, levels = design)

# Apply the contrast to the model fit
fit2 <- contrasts.fit(fit, contrast_matrix)

# Step 6: Apply eBayes to borrow strength across genes
fit2 <- eBayes(fit2)

# Step 7: Find the gene with the lowest p-value
top_genes <- topTable(fit2, number = 1, sort.by = "P")

# Step 8: Extract the absolute value of the log fold change (logFC)
abs_logFC <- abs(top_genes$logFC)

# Output the result
cat("Absolute value of the logFC for the gene with the lowest p-value:", abs_logFC, "\n")


## Question 7 

# Step 1: Extract the topTable with all genes (to access adj.P.value)
all_genes <- topTable(fit2, number = Inf)  # Extract all genes

# Step 2: Count how many genes have an adj.P.value < 0.05
num_diff_expressed_genes <- sum(all_genes$adj.P.Val < 0.05)

# Output the result
cat("Number of differentially expressed genes with adj.P.value < 0.05:", num_diff_expressed_genes, "\n")


## Question 8

# Step 1: Install and load necessary packages
if (!requireNamespace("minfi", quietly = TRUE)) {
  BiocManager::install("minfi")
}
if (!requireNamespace("minfiData", quietly = TRUE)) {
  BiocManager::install("minfiData")
}

# Load the necessary libraries
library(minfi)
library(minfiData)

# Step 2: Load the RGsetEx dataset
data(RGsetEx)

# Step 3: Preprocess the data using the preprocessFunnorm function
preprocessed_data <- preprocessFunnorm(RGsetEx)

# Step 4: Extract Beta values (percent methylation)
beta_values <- getBeta(preprocessed_data)

# Step 5: Identify OpenSea loci (assumed to be in feature data or annotation)
# Get annotation and filter OpenSea CpGs (use 'getAnnotation()' to check CpG islands, shores, etc.)
annotation_data <- getAnnotation(preprocessed_data)

# Filter for OpenSea loci (annotation$Relation_to_Island == "OpenSea")
open_sea_indices <- which(annotation_data$Relation_to_Island == "OpenSea")
beta_open_sea <- beta_values[open_sea_indices, ]

# Step 6: Compute the average Beta values for each sample across OpenSea CpGs
avg_beta_per_sample <- colMeans(beta_open_sea, na.rm = TRUE)

# Step 7: Define groups (3 normals and 3 cancer samples based on the dataset structure)
# Assuming the first 3 samples are normal and the last 3 samples are cancer
normal_samples <- avg_beta_per_sample[1:3]
cancer_samples <- avg_beta_per_sample[4:6]

# Step 8: Calculate the mean difference in Beta values between normal and cancer samples
mean_diff_beta <- mean(normal_samples) - mean(cancer_samples)

# Output the result
cat("Mean difference in Beta values between normal and cancer samples (OpenSea CpGs):", mean_diff_beta, "\n")


## Question 9 

# Load the required libraries
library(AnnotationHub)
library(minfi)
library(GenomicRanges)

# Step 2: Access AnnotationHub to get DNase hypersensitive sites for the Caco2 cell line
ah <- AnnotationHub()

# Query for DNase hypersensitive sites for the Caco2 cell line from ENCODE
query_results <- query(ah, c("Caco2", "DNase", "narrowPeak", "ENCODE"))
dnase_data <- query_results[["AH22442"]]

# Step 3: Load the 450k dataset (assuming it's the same RGsetEx dataset used in Question 8)
data(RGsetEx)

# Get the annotation for the 450k array
annotation_data <- getAnnotation(RGsetEx)

# Step 4: Create a GRanges object for the CpGs on the 450k array
cpg_granges <- GRanges(seqnames = annotation_data$chr,
                       ranges = IRanges(start = annotation_data$pos, end = annotation_data$pos))

# Step 5: Create a GRanges object for the DNase hypersensitive sites
dnase_granges <- GRanges(seqnames = seqnames(dnase_data),
                         ranges = IRanges(start = start(dnase_data), end = end(dnase_data)))

# Step 6: Find overlaps between DNase hypersensitive sites and CpGs on the 450k array
overlaps <- findOverlaps(dnase_granges, cpg_granges)

# Step 7: Count how many unique DNase sites contain one or more CpGs
dnase_with_cpgs <- unique(queryHits(overlaps))
num_dnase_with_cpgs <- length(dnase_with_cpgs)

# Output the result
cat("Number of DNase hypersensitive sites containing one or more CpGs on the 450k array:", num_dnase_with_cpgs, "\n")


## Question 10

# Step 1: Install and load the necessary packages
if (!requireNamespace("zebrafishRNASeq", quietly = TRUE)) {
  BiocManager::install("zebrafishRNASeq")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

# Load the required libraries
library(zebrafishRNASeq)
library(DESeq2)

# Step 2: Load the zfGenes dataset
data(zfGenes)

# Step 3: Exclude rows with row names starting with "ERCC" (spike-ins)
zfGenes_filtered <- zfGenes[!grepl("^ERCC", rownames(zfGenes)), ]

# Step 4: Set up sample information for DESeq2 (control and treatment groups)
# Assuming the column names in zfGenes correspond to sample names
sample_info <- data.frame(
  condition = factor(rep(c("control", "treatment"), each = 3))
)
rownames(sample_info) <- colnames(zfGenes_filtered)

# Step 5: Create a DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = zfGenes_filtered, colData = sample_info, design = ~ condition)

# Step 6: Run the DESeq2 differential expression analysis
dds <- DESeq(dds)

# Step 7: Get the results with p-adjusted values
results_DESeq2 <- results(dds)

# Step 8: Filter for differentially expressed features (padj <= 0.05)
diff_expressed <- results_DESeq2[!is.na(results_DESeq2$padj) & results_DESeq2$padj <= 0.05, ]

# Step 9: Count how many features are differentially expressed
num_diff_expressed <- nrow(diff_expressed)

# Output the result
cat("Number of differentially expressed features (padj <= 0.05):", num_diff_expressed, "\n")
