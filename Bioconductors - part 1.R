library(Biostrings)

dna2 <- DNAString("ACGT-N")
letterFrequency(dna2, "GC")

BiocManager::install("rtracklayer")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the downloaded package
install.packages("BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz", repos = NULL, type = "source")
# Load necessary libraries
library(Biostrings)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Load the hg19 genome data
hg19_genome <- BSgenome.Hsapiens.UCSC.hg19

# Extract chromosome 22 sequence from hg19
chr22_seq <- hg19_genome$chr22

# Convert DNA sequence to a character string
chr22_seq_string <- as.character(chr22_seq)

# Remove "N" bases using gsub
chr22_clean_string <- gsub("N", "", chr22_seq_string)

# Convert back to DNAString object
chr22_clean <- DNAString(chr22_clean_string)


# Calculate the frequency of G and C bases
gc_content <- letterFrequency(chr22_clean, "GC", as.prob = TRUE)

# Print the GC content as a percentage
cat("GC content of chr22 in hg19:", gc_content * 100, "%\n")



## Question 2 
# Install and load necessary packages
# Install the AnnotationHub package if not already installed
if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
  BiocManager::install("AnnotationHub")
}

# Load the AnnotationHub package
library(AnnotationHub)

# Initialize AnnotationHub
ah <- AnnotationHub()

# Query for H3K27me3 narrowPeak data from E003 (H1 stem cell line)
query_results <- query(ah, c("E003", "H3K27me3", "narrowPeak"))

# Step 1: Retrieve the narrowPeak data for H3K27me3 in E003 (H1 stem cell line)
narrowPeak_data <- query_results[["AH29892"]]

# Step 2: Filter the narrowPeak data for chromosome 22
chr22_peaks <- narrowPeak_data[seqnames(narrowPeak_data) == "chr22"]

# Filter narrowPeak data for chr22
chr22_peaks <- narrowPeak_data[seqnames(narrowPeak_data) == "chr22"]



# Function to calculate GC content as a percentage
calculate_gc_content <- function(seq) {
  gc_content <- letterFrequency(seq, letters = "GC", as.prob = TRUE)
  return(gc_content * 100)  # Return GC content as percentage
}

# Calculate GC content for each peak on chr22
gc_contents <- sapply(seq_along(chr22_peaks), function(i) {
  start_pos <- start(chr22_peaks[i])
  end_pos <- end(chr22_peaks[i])
  peak_seq <- subseq(chr22_seq, start = start_pos, end = end_pos)
  calculate_gc_content(peak_seq)
})

# Compute the average GC content across all peaks
mean_gc_content <- mean(gc_contents)

# Output the mean GC content
cat("Mean GC content of H3K27me3 narrowPeak regions on chr22:", mean_gc_content / 100, "\n")

## Question 3 

# Example URL for H3K27me3 narrowPeak file for E003
peak_file_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E003-H3K27me3.narrowPeak.gz"

# Download the file
download.file(peak_file_url, destfile = "E003-H3K27me3.narrowPeak.gz")

# Load necessary libraries
library(GenomicRanges)
library(rtracklayer)

# Read the narrowPeak file
peaks <- import("E003-H3K27me3.narrowPeak.gz", format = "narrowPeak")

# Filter for chr22 peaks
chr22_peaks <- peaks[seqnames(peaks) == "chr22"]

# Install and load the hg19 genome package if needed
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

# Load the hg19 genome
library(BSgenome.Hsapiens.UCSC.hg19)
hg19_genome <- BSgenome.Hsapiens.UCSC.hg19

# Extract the sequence of chr22
chr22_seq <- hg19_genome$chr22


# Function to calculate GC content as a percentage
calculate_gc_content <- function(start, end) {
  seq <- subseq(chr22_seq, start=start, end=end)
  gc_content <- letterFrequency(seq, "GC", as.prob = TRUE)
  return(gc_content * 100)  # Return GC content as a percentage
}

# Calculate GC content for each peak on chr22
gc_contents <- mapply(calculate_gc_content, start(chr22_peaks), end(chr22_peaks))


# Extract the signalValue column
signal_values <- mcols(chr22_peaks)$signalValue

# Calculate the correlation between GC content and signalValue
gc_signal_correlation <- cor(gc_contents, signal_values, use = "complete.obs")

# Output the correlation
cat("Correlation between GC content and signalValue on chr22:", gc_signal_correlation, "\n")


## Question 4 

# Install and load AnnotationHub
if (!requireNamespace("AnnotationHub", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub")

library(AnnotationHub)
library(rtracklayer)

# Create an AnnotationHub object to search for data
ah <- AnnotationHub()

# Query for H3K27me3 data from E003 (H1 cell line) in the hub
query_result <- query(ah, c("E003", "H3K27me3", "fc.signal"))

# Display the query result to find the appropriate dataset
query_result


# Load the fc.signal data using its AnnotationHub ID (assuming the ID is found in the query result)
fc_signal_data <- query_result[["AH32033"]]

# Print a summary of the loaded data
summary(fc_signal_data)

# Check the type of the fc_signal_data object
class(fc_signal_data)

# Import necessary libraries
library(GenomicRanges)
library(rtracklayer)

# Import fc.signal data for chr22
chr22_range <- GRanges("chr22", IRanges(start = 1, end = 51304566))  # Length of chr22 in hg19
fc_signal_chr22 <- import(fc_signal_data, which = chr22_range)

# Define a function to compute the average fc.signal for each peak
compute_avg_fc_signal <- function(peak, signal_data) {
  start_pos <- start(peak)
  end_pos <- end(peak)
  
  # Extract the signal data (fc.signal) for this region
  region_signal <- subsetByOverlaps(signal_data, peak)$score
  
  # Return the mean fc.signal for the region
  return(mean(region_signal, na.rm = TRUE))
}

# Compute the average fc.signal for each peak
avg_fc_signal <- sapply(seq_along(chr22_peaks), function(i) {
  compute_avg_fc_signal(chr22_peaks[i], fc_signal_chr22)
})

# Extract the signalValue column from the narrowPeak regions
signal_values <- mcols(chr22_peaks)$signalValue

# Calculate the correlation between signalValue and average fc.signal
correlation <- cor(signal_values, avg_fc_signal, use = "complete.obs")

# Output the correlation
cat("Correlation between signalValue and average fc.signal:", correlation, "\n")


## Question 5 

library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Question 5: Count how many bases on chr22 have an fc.signal >= 1

# Import fc.signal data for chr22
chr22_range <- GRanges("chr22", IRanges(start = 1, end = 51304566))  # Length of chr22 in hg19
fc_signal_chr22 <- import(fc_signal_data, which = chr22_range)



#  How many bases on chr22 have an fc.signal >= 1?
# Find the number of bases where fc.signal >= 1
fc_signal_values <- as.numeric(fc_signal_chr22$score)
bases_with_fc_signal_gte_1 <- sum(fc_signal_values >= 1)

# Print the result
cat("Number of bases on chr22 with fc.signal >= 1:", bases_with_fc_signal_gte_1, "\n")
                 

## Question 6 

# Load necessary libraries
library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Initialize AnnotationHub
ah <- AnnotationHub()

# Query for H3K27me3 fc.signal for E003 (H1 stem cell line)
query_e003 <- query(ah, c("E003", "H3K27me3", "fc.signal"))

# Query for H3K27me3 fc.signal for E055 (foreskin fibroblast)
query_e055 <- query(ah, c("E055", "H3K27me3", "fc.signal"))

# Load fc.signal data for E003 (H1 stem cell line)
fc_signal_e003 <- query_e003[["AH32033"]]

# Load fc.signal data for E055 (foreskin fibroblast)
fc_signal_e055 <- query_e055[["AH32470"]]

# Filter fc.signal data for chr22
# Import fc.signal data for chr22
chr22_range <- GRanges("chr22", IRanges(start = 1, end = 51304566))  # Length of chr22 in hg19
#fc_signal_chr22 <- import(fc_signal_data, which = chr22_range)

fc_signal_e003_chr22 <- import(fc_signal_e003, which = chr22_range)
fc_signal_e055_chr22 <- import(fc_signal_e055, which = chr22_range)



# Filter for regions where E003 signal <= 0.5
ir_e003_filtered <- fc_signal_e003_chr22[mcols(fc_signal_e003_chr22)$score <= 0.5]

# Filter for regions where E055 signal >= 2
ir_e055_filtered <- fc_signal_e055_chr22[mcols(fc_signal_e055_chr22)$score >= 2]

# Find intersecting regions where the conditions for both E003 and E055 are satisfied
overlap_regions <- intersect(ir_e003_filtered, ir_e055_filtered)

# Output the overlapping regions
overlap_regions


## Question 7 


# Load necessary libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(Biostrings)

# Step 1: Load hg19 genome and get chr22 sequence
hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
chr22_seq <- hg19_genome$chr22



# Step 2: Load CpG island data from UCSC for hg19 and filter for chr22
cpg_islands_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz"
download.file(cpg_islands_url, destfile = "cpgIslandExt.txt.gz")
cpg_islands <- read.table(gzfile("cpgIslandExt.txt.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(cpg_islands) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpG", "perGC", "obsExp")

cpg_islands_chr22 <- subset(cpg_islands, chrom == "chr22")

library(GenomicRanges)
cpg_islands_granges <- GRanges(seqnames = cpg_islands_chr22$chrom,
                               ranges = IRanges(start = cpg_islands_chr22$chromStart,
                                                end = cpg_islands_chr22$chromEnd),
                               strand = "*")


# Step 3: Define function to calculate observed-to-expected CpG ratio
calculate_cpg_ratio <- function(seq) {
  # Use countPattern for single sequence
  observed_cpg <- countPattern("CG", seq)
  
  # Get frequencies of C and G
  c_count <- letterFrequency(seq, "C", as.prob = TRUE)
  g_count <- letterFrequency(seq, "G", as.prob = TRUE)
  
  # Calculate expected CpG
  expected_cpg <- (c_count * g_count) * nchar(seq)
  
  # Return the observed-to-expected CpG ratio
  if (expected_cpg > 0) {
    return(observed_cpg / expected_cpg)
  } else {
    return(NA)
  }
}

# Step 4: Calculate CpG ratios for each CpG island on chr22 using the GRanges object
cpg_ratios <- sapply(seq_along(cpg_islands_granges), function(i) {
  # Get start and end positions of the i-th CpG island
  start_pos <- start(cpg_islands_granges[i])
  end_pos <- end(cpg_islands_granges[i])
  
  # Extract the corresponding sequence from chr22
  region_seq <- subseq(chr22_seq, start = start_pos, end = end_pos)
  
  # Calculate CpG ratio for the sequence
  calculate_cpg_ratio(region_seq)
})

# Step 5: Calculate average observed-to-expected CpG ratio
average_cpg_ratio <- mean(cpg_ratios, na.rm = TRUE)
cat("Average observed-to-expected CpG ratio:", average_cpg_ratio, "\n")

## Question 8 

# Load necessary libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

# Step 1: Load hg19 genome and get chr22 sequence
#hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
chr22_seq <- hg19_genome$chr22

# Step 2: Search for "TATAAA" on the forward strand and "TTTATA" on the reverse strand
tata_count_forward <- countPattern("TATAAA", chr22_seq)
tata_count_reverse <- countPattern("TTTATA", chr22_seq)

# Step 3: Total TATA boxes on both strands
total_tata_boxes <- tata_count_forward + tata_count_reverse

# Output the total number of TATA boxes
cat("Total number of TATA boxes on chr22:", total_tata_boxes, "\n")


## Question 9 

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Biostrings)

# Step 1: Load transcript and CDS information for chr22
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Load all transcripts
transcripts_all <- transcripts(txdb, columns = c("tx_id", "tx_name"))

# Filter for chr22
transcripts_chr22 <- subset(transcripts_all, seqnames(transcripts_all) == "chr22")

# Load all CDS regions
cds_all <- cds(txdb, columns = c("cds_id", "tx_id"))

# Filter for chr22
cds_chr22 <- subset(cds_all, seqnames(cds_all) == "chr22")

# Step 2: Define promoter regions (900 bp upstream, 100 bp downstream)
promoters_chr22 <- promoters(transcripts_chr22, upstream = 900, downstream = 100)

mcols(promoters_chr22)$tx_id <- mcols(transcripts_chr22)$tx_id

# Step 3: Filter promoters for transcripts with coding sequences
cds_tx_ids <- unique(mcols(cds_chr22)$tx_id)
# Unlist cds_tx_ids to convert it into a regular integer vector
cds_tx_ids_vector <- unlist(cds_tx_ids)

# Filter promoters for transcripts that have coding sequences
promoters_with_cds <- promoters_chr22[mcols(promoters_chr22)$tx_id %in% cds_tx_ids_vector]


# Step 4: Load chr22 sequence
#hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
chr22_seq <- hg19_genome$chr22

# Step 5: Extract promoter sequences from chr22 genome
#promoter_sequences <- getSeq(hg19_genome, promoters_with_cds)

promoter_sequences <- getSeq(hg19_genome, promoters_with_cds)

# Step 4: Define function to search for TATA boxes based on strand and use lapply
contains_tata <- lapply(seq_along(promoter_sequences), function(i) {
  seq <- promoter_sequences[[i]]
  # Convert the strand to a character
  strand <- as.character(promoter_strands[i])
  
  if (strand == "+") {
    return(countPattern("TATAAA", seq, fixed = TRUE))
  } else if (strand == "-") {
    return(countPattern("TTTATA", seq, fixed = TRUE))
  } else {
    return(0)
  }
})








# Filter promoters that contain TATA boxes on the correct strand
promoters_with_tata <- promoters_with_cds[contains_tata > 0]

# Step 6: Count the number of promoters with TATA boxes
num_promoters_with_tata <- length(promoters_with_tata)
cat("Number of promoters on chr22 with TATA boxes on the same strand as the transcript:", num_promoters_with_tata, "\n")

## Question 10 

# Load required libraries
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

# Step 1: Load transcript and CDS information for chr22
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Load all transcripts
transcripts_all <- transcripts(txdb, columns = c("tx_id", "tx_name"))

# Filter for chr22
transcripts_chr22 <- subset(transcripts_all, seqnames(transcripts_all) == "chr22")

# Load all CDS regions
cds_all <- cds(txdb, columns = c("cds_id", "tx_id"))

# Filter for chr22
cds_chr22 <- subset(cds_all, seqnames(cds_all) == "chr22")

# Step 2: Define promoter regions (900 bp upstream, 100 bp downstream)
promoters_chr22 <- promoters(transcripts_chr22, upstream = 900, downstream = 100)

mcols(promoters_chr22)$tx_id <- mcols(transcripts_chr22)$tx_id

# Step 3: Filter promoters for transcripts with coding sequences
cds_tx_ids <- unique(mcols(cds_chr22)$tx_id)
# Unlist cds_tx_ids to convert it into a regular integer vector
cds_tx_ids_vector <- unlist(cds_tx_ids)

# Filter promoters for transcripts that have coding sequences
promoters_with_cds <- promoters_chr22[mcols(promoters_chr22)$tx_id %in% cds_tx_ids_vector]


# Step 4: Find overlapping promoter regions, excluding self-overlaps
overlaps <- findOverlaps(promoters_with_cds)

# Filter out self-overlaps (where a promoter is overlapping itself)
valid_overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]

# Extract only the actual overlapping regions using pintersect
overlapping_regions <- pintersect(promoters_with_cds[queryHits(valid_overlaps)],
                                  promoters_with_cds[subjectHits(valid_overlaps)])

# Merge overlapping regions using reduce to avoid double-counting
reduced_overlapping_promoters <- reduce(overlapping_regions)

# Step 5: Count the number of unique bases in overlapping regions
total_bases_overlapping <- sum(width(reduced_overlapping_promoters))

# Output the result
cat("Number of bases on chr22 that are part of more than one promoter:", total_bases_overlapping, "\n")
