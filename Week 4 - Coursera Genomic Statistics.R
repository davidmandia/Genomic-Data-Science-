## Question 1 

BiocManager::install("goseq")

# Load the goseq package
library(goseq)

# Get the list of supported genomes in goseq
supported_genomes <- supportedGenomes()

# Print the supported Mouse genome builds
supported_genomes[grep("mm", supported_genomes)]

## Question 2 

library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
edata = log2(edata+1)


# Load the Bottomly data with the following code and perform a differential expression analysis using 
# limma with only the strain variable as an outcome. 
# How many genes are differentially expressed at the 5% FDR level using Benjamini-Hochberg correction?
# What is the gene identifier of the first gene differentially expressed at this level
# (just in order, not the smallest FDR)
# ? (hint: the featureNames function may be useful)

# Load necessary libraries
library(limma)

# Create a design matrix using the 'strain' variable
design <- model.matrix(~ pdata_bot$strain)


# Fit the linear model using lmFit
fit <- lmFit(edata, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract the top table of results
limma_res <- topTable(fit, coef = 2, number = nrow(edata))  # coef = 2 refers to the strain variable

# Apply Benjamini-Hochberg correction to the p-values
limma_res$adj.P.Val <- p.adjust(limma_res$P.Value, method = "BH")

# Filter genes that are significant at 5% FDR level
significant_genes <- limma_res[limma_res$adj.P.Val < 0.05, ]

# Number of significant genes
num_significant_genes <- nrow(significant_genes)
cat("Number of differentially expressed genes at 5% FDR level:", num_significant_genes, "\n")

# Get the gene identifier of the first differentially expressed gene
first_gene <- rownames(significant_genes)[1]

# Display the gene identifier
cat("The gene identifier of the first differentially expressed gene at 5% FDR:", first_gene, "\n")

## Question 3 

# Create a binary vector where 1 = differentially expressed, 0 = not
gene_list <- as.integer(rownames(edata) %in% rownames(significant_genes))
names(gene_list) <- rownames(edata)


# Perform bias estimation (nullp) using Ensembl gene IDs (ensGene) for the mm9 genome
pwf <- nullp(gene_list, "mm9", "ensGene")

## required for goseq
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
# Perform GO analysis
go_results <- goseq(pwf, "mm9", "ensGene")

# Display the top over-represented categories
head(go_results[order(go_results$over_represented_pvalue), ])


# Find the top GO category based on the over-represented p-value
top_category <- go_results[which.min(go_results$over_represented_pvalue), ]

# Print the top category
cat("Top over-represented GO category:", top_category$category, "\n")
cat("Description:", top_category$term, "\n")
cat("P-value:", top_category$over_represented_pvalue, "\n")


## Question 5 

# Load the Bottomly data with the following code and perform a differential expression analysis using 
# limma and treating strain as the outcome but adjusting for lane as a factor.
# Then find genes significant at the 5% FDR rate using the Benjamini Hochberg correction
#and perform the gene set analysis with goseq following the protocol from the first 4 questions.
#How many of the top 10 overrepresented categories are the same for the adjusted and unadjusted analysis?


library(Biobase)
library(limma)

# Load the Bottomly dataset
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)

# Extract data
bot = bottomly.eset
pdata_bot = pData(bot)
edata = exprs(bot)

# Filter out lowly expressed genes
edata = edata[rowMeans(edata) > 5, ]

# Create a design matrix including strain and lane
design_adjusted <- model.matrix(~ pdata_bot$strain + pdata_bot$lane)


# Fit the linear model with adjustment for lane
fit_adjusted <- lmFit(edata, design_adjusted)
fit_adjusted <- eBayes(fit_adjusted)

# Get the top genes with strain as the outcome, adjusting for lane
limma_res_adjusted <- topTable(fit_adjusted, coef = 2, number = nrow(edata))  # coef = 2 corresponds to strain


# Apply Benjamini-Hochberg correction to the p-values
limma_res_adjusted$adj.P.Val <- p.adjust(limma_res_adjusted$P.Value, method = "BH")

# Filter significant genes at 5% FDR
significant_genes_adjusted <- limma_res_adjusted[limma_res_adjusted$adj.P.Val < 0.05, ]


# Create a binary vector of differentially expressed genes
gene_list_adjusted <- as.integer(rownames(edata) %in% rownames(significant_genes_adjusted))
names(gene_list_adjusted) <- rownames(edata)

# Perform bias estimation (nullp) for the mm9 genome
pwf_adjusted <- nullp(gene_list_adjusted, "mm9", "ensGene")

# Perform GO analysis
go_results_adjusted <- goseq(pwf_adjusted, "mm9", "ensGene")

# Get the top 10 over-represented GO categories (sorted by p-value)
top10_adjusted <- go_results_adjusted[order(go_results_adjusted$over_represented_pvalue), ][1:10, ]


## Not adjusted for strain 
# Create a design matrix with only strain
design_unadjusted <- model.matrix(~ pdata_bot$strain)

# Fit the linear model without adjusting for lane
fit_unadjusted <- lmFit(edata, design_unadjusted)
fit_unadjusted <- eBayes(fit_unadjusted)

# Get the top genes
limma_res_unadjusted <- topTable(fit_unadjusted, coef = 2, number = nrow(edata))

# Apply Benjamini-Hochberg correction
limma_res_unadjusted$adj.P.Val <- p.adjust(limma_res_unadjusted$P.Value, method = "BH")

# Filter significant genes at 5% FDR
significant_genes_unadjusted <- limma_res_unadjusted[limma_res_unadjusted$adj.P.Val < 0.05, ]

# Create a binary vector for significant genes (unadjusted)
gene_list_unadjusted <- as.integer(rownames(edata) %in% rownames(significant_genes_unadjusted))
names(gene_list_unadjusted) <- rownames(edata)

# Perform bias estimation for unadjusted analysis
pwf_unadjusted <- nullp(gene_list_unadjusted, "mm9", "ensGene")

# Perform GO analysis for unadjusted results
go_results_unadjusted <- goseq(pwf_unadjusted, "mm9", "ensGene")

# Get the top 10 over-represented GO categories for unadjusted analysis
top10_unadjusted <- go_results_unadjusted[order(go_results_unadjusted$over_represented_pvalue), ][1:10, ]

# Find how many of the top 10 categories are the same between the adjusted and unadjusted analyses
common_categories <- intersect(top10_adjusted$category, top10_unadjusted$category)

# Count the number of common categories
num_common_categories <- length(common_categories)

# Output the result
cat("Number of common top 10 over-represented categories:", num_common_categories, "\n")

