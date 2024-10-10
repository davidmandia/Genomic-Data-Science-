BiocManager::install("snpStats")
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc


#fit a linear model and a logistic regression model to the data for the 3rd SNP.
#What are the coefficients for the SNP variable?
#How are they interpreted? (Hint: Don't forget to recode the 0 values to NA for the SNP data)

# Convert snpdata to a data frame if it's not already
snpdata_df <- as.data.frame(snpdata)

# Recode 0 values to NA using apply
snpdata_df <- apply(snpdata_df, 2, function(x) { x[x == 0] <- NA; return(x) })

# Fit the models using the updated data frame
linear_model <- lm(status ~ snpdata_df[, 3], na.action = na.exclude)
logistic_model <- glm(status ~ snpdata_df[, 3], family = binomial, na.action = na.exclude)

# Get coefficients
# For the linear model
linear_coefficients_third_snp <- summary(linear_model)$coefficients[3, ]
print(linear_coefficients_third_snp)

# For the logistic model
logistic_coefficients_third_snp <- summary(logistic_model)$coefficients[3, ]
print(logistic_coefficients_third_snp)


## Question 3


library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc


#Fit a logistic regression model on a recessive
#(need 2 copies of minor allele to confer risk)
#and additive scale for the 10th SNP.
#Make a table of the fitted values versus the case/control status.
#Does one model fit better than the other?

dim(snpdata)

snp_10 = as.numeric(snpdata[,10])
dim(snp_10)


table(snp_10)

snp_10[snp_10 == 0] = NA

# estimate is the change in log odds per additional allele 
glm1  = glm(status ~ snp_10, family = "binomial")

tidy(glm1)


snp_10_dom = (snp_10 == 2)

glm2 = glm(status ~ snp_10_dom, family = "binomial")

tidy(glm2)
