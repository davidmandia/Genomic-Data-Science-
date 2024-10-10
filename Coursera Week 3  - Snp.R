install.packages(c("devtools","broom","MASS"))


install.packages("BiocManager")
BiocManager::install(c("Biobase", "snpStats", "DESeq2"))
install.packages("broom")



## Question 1
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc


xxmat <- xxt(sub.10, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]


## Fit a linear model and a logistic regression model to the data for the 3rd SNP.
## What are the coefficients for the SNP variable? How are they interpreted? 
## (Hint: Don't forget to recode the 0 values to NA for the SNP data)

snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA

## Logistic  - coeff is the decrease in log odds ration with each additional copy
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3)

## Linear  -  coeff is the decrease in proba for each additional copy
fit_lm_3 <- lm(status ~ snp3)
tidy(fit_lm_3)


## Question 3

use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc


## Fit a logistic regression model on a recessive (need 2 copies of minor allele to confer risk) 
## and additive scale for the 10th SNP. Make a table of the fitted values versus the case/control status.
## Does one model fit better than the other?

snp_10 <- as.numeric(snp_10)
# Recode for the recessive model:
# 0: Homozygous major or heterozygous
# 1: Homozygous minor
snp_10_rec <- ifelse(snp_10 == 2, 1, 0)


# Fit the logistic regression model under the recessive model
fit_logit_rec <- glm(status ~ snp_10_rec, family = binomial(link = "logit"))

# View summary of the model
summary(fit_logit_rec)


# Fit the logistic regression model under the additive model
fit_logit_add <- glm(status ~ snp_10, family = binomial(link = "logit"))

# View summary of the model
summary(fit_logit_add)


# Predicted probabilities for the recessive model
fitted_rec <- fitted(fit_logit_rec)

# Predicted probabilities for the additive model
fitted_add <- fitted(fit_logit_add)

# Create tables to compare fitted values to case/control status
table_rec <- table(round(fitted_rec), status)
table_add <- table(round(fitted_add), status)

# Print the tables
table_rec
table_add


## Question 4 

library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc


## Fit an additive logistic regression model to each SNP.
## What is the average effect size? What is the max? What is the minimum?


#snpdata = as.numeric(snpdata)
#snpdata[snpdata==0] = NA

# Initialize a vector to store the effect sizes (coefficients)
effect_sizes <- numeric(ncol(snpdata))





# Loop over each SNP in the snpdata matrix (each column is a SNP)
for (i in 1:ncol(snpdata)) {
  # Extract the ith SNP genotype data
  snp_i <- as.numeric(snpdata[, i])
  
  # Fit a logistic regression model assuming an additive effect
  fit <- glm(status ~ snp_i, family = binomial(link = "logit"))
  
  # Extract the coefficient for the SNP (2nd element of coefficients: snp_i)
  effect_sizes[i] <- coef(fit)[2]
}

# Display the coefficients (effect sizes)
effect_sizes


# Calculate the average effect size
average_effect_size <- mean(effect_sizes)

# Calculate the maximum effect size
max_effect_size <- max(effect_sizes)

# Calculate the minimum effect size
min_effect_size <- min(effect_sizes)

# Display the results
average_effect_size
max_effect_size
min_effect_size



## Question 5 

library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc


## Fit an additive logistic regression model to each SNP and square the coefficients.
## What is the correlation with the results from using  snp.rhs.tests and chi.squared?
# Why does this make sense?

glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)
slotNames(glm_all)


# Summary of the test results (chi-squared statistics and p-values)
summary(glm_all)

# Extract p-values
p_values <- p.value(glm_all)

qq.chisq(chi.squared(glm_all),df=1)

# Initialize a vector to store the squared coefficients
squared_coeffs <- numeric(ncol(snpdata))

# Loop over each SNP in the snpdata matrix
for (i in 1:ncol(snpdata)) {
  # Extract the ith SNP genotype data
  snp_i <- as.numeric(snpdata[, i])
  
  # Fit an additive logistic regression model
  fit <- glm(status ~ snp_i, family = binomial(link = "logit"))
  
  # Extract the coefficient for the SNP and square it (second coefficient: snp_i)
  squared_coeffs[i] <- coef(fit)[2]^2
}


# Extract the chi-squared statistics from the snp.rhs.tests result
chisq_values <- chi.squared(glm_all)

# Calculate the correlation between the squared coefficients and chi-squared statistics
correlation <- cor(squared_coeffs, chisq_values)

# Print the correlation
correlation


## Question 6 

library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

## Do the log2(data + 1) transform and fit calculate F-statistics
## for the difference between studies/populations using genefilter:rowFtests and using genefilter:rowttests.
## Do you get the same statistic? Do you get the same p-value?

edata = log2(as.matrix(edata) + 1)


# Define the group variable (assuming 'group' in pdata represents population/study)
group <- pdata$study
install.packages("devtools")

BiocManager::install(c("Biobase", "limma", "genefilter"))
BiocManager::install("genefilter")
library(genefilter)


## T-test 
tstats_obj = rowttests(edata,group)

names(tstats_obj)
hist(tstats_obj$statistic,col=2)


## Row F-test 

fstats_obj = rowFtests(edata, group)
names(fstats_obj)

hist(tstats_obj$statistic,col=2)


# Question 7

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)


# First test for differences between the studies using the DESeq2 package using the DESeq function.
# Then do the log2(data + 1) transform and do the test for differences between studies using the limma  package
# and the lmFit, ebayes and topTable functions.
# What is the correlation in the statistics between the two analyses? 
# Are there more differences for the large statistics or the small statistics (hint: Make an MA-plot).

BiocManager::install("DESeq2")

# Load the DESeq2 package
library(DESeq2)


# Create a DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = edata,
                              colData = pdata,
                              design = ~ study)  # Assuming "study" indicates the different groups

# Run the DESeq analysis
dds <- DESeq(dds)

# Extract the results
res_deseq2 <- results(dds)

# Get the log2 fold changes and p-values from DESeq2
deseq2_stat <- res_deseq2$stat
deseq2_pval <- res_deseq2$pvalue


# Install limma if needed
BiocManager::install("limma")

# Load limma
library(limma)


# Perform log2 transformation on expression data
log_edata <- log2(edata + 1)

# Create a design matrix for the linear model
design <- model.matrix(~ 0 + pdata$study)  # Assuming "study" indicates the groups

# Fit the linear model to the log2-transformed data
fit <- lmFit(log_edata, design)

# Apply empirical Bayes smoothing to the standard errors
fit <- eBayes(fit)

# Extract the top table (differential expression results)
limma_res <- topTable(fit, number = nrow(log_edata))

# Get the statistics and p-values from limma
limma_stat <- limma_res$
limma_pval <- limma_res$P.Value

length(deseq2_stat)
length(limma_stat)



# Calculate the correlation between DESeq2 and limma statistics
correlation <- cor(deseq2_stat, limma_stat, use = "complete.obs")

# Print the correlation result
print(correlation)


# Check for NA values in log2FoldChange and baseMean
na_log2FoldChange <- sum(is.na(res_deseq2$log2FoldChange))
na_baseMean <- sum(is.na(res_deseq2$baseMean))

# Print the results
cat("Number of NA values in log2FoldChange:", na_log2FoldChange, "\n")
cat("Number of NA values in baseMean:", na_baseMean, "\n")


# Create an MA-plot for DESeq2
plotMA(res_deseq2, main = "MA-plot DESeq2", ylim = c(-2, 2))


# Create an MA-plot for limma
plotMA(fit, coef = 1, main = "MA-plot limma", ylim = c(-2, 2))


## Question 8 

# Apply Benjamini-Hochberg correction to DESeq2 p-values
bh_adjusted_pvals_deseq2 <- p.adjust(res_deseq2$pvalue, method = "BH")

# Count how many p-values are significant at FDR 0.05
significant_deseq2 <- sum(bh_adjusted_pvals_deseq2 < 0.05)

# Print the number of significant results
cat("Number of statistically significant results in DESeq2 at FDR 0.05:", significant_deseq2, "\n")


# Apply Benjamini-Hochberg correction to limma p-values
bh_adjusted_pvals_limma <- p.adjust(limma_res$P.Value, method = "BH")

# Count how many p-values are significant at FDR 0.05
significant_limma <- sum(bh_adjusted_pvals_limma < 0.05)

# Print the number of significant results
cat("Number of statistically significant results in limma at FDR 0.05:", significant_limma, "\n")


