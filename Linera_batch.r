install.packages('BiocManager')

BiocManager::install("Biobase")
BiocManager::install("GenomicRanges")

install.packages("devtools")
library(devtools)

library(Biobase)
library(GenomicRanges)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# 1st principal component in the data set if you

# Do no transformations?
# log2(data + 1) transform?
# log2(data + 1) transform and subtract row means?

svd_no = svd(edata)

svd_no$d^2/sum(svd_no$d^2)*100

#b
svd_log = svd(log2(edata+1))

svd_log$d^2/sum(svd_log$d^2)*100

#c 
edata = log2(edata +1)


edata_centered = edata -rowMeans(edata)

svd_row = svd(edata_centered)


rowmeans = rowMeans(edata)

svd_row$d^2/sum(svd_row$d^2)*100


## question 2

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# What is the correlation between the first singular vector and the sample clustering indicator?
log_data = log2(edata + 1)

edata_centered = edata -rowMeans(edata)
mean_centered_data = sweep(log_data, 1, rowMeans(log_data))

set.seed(333)


kmeans1 = kmeans(mean_centered_data, centers = 2)

svd1 = svd(mean_centered_data)


cluster_indicator <- kmeans1$cluster

first_singular_vector <- svd1$u[, 1]

correlation <- cor(first_singular_vector, as.numeric(cluster_indicator))
print(correlation)

head(mean_centered_data)
table(cluster_indicator)
str(svd1)
correlation <- cor(first_singular_vector, as.numeric(cluster_indicator))
print(correlation)
sum(is.na(mean_centered_data))  # Check for NAs
is.numeric(mean_centered_data)   # Check if it's numeric



### Question 3

library(broom)
library(devtools)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

table(pdata_bm$num.tech.reps)

first_gene = edata[1,]

lm1 = lm(first_gene ~ pdata_bm$num.tech.reps)

tidy(lm1)



## Question 7 

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

### Question 7
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

#Fit many regression models to the expression
#data where age age is the outcome variable
#using the lmFit lmFit function from the limma
#limma package (hint: you may have to subset
#the expression data to the samples without 
#missing values of age to get the model to fit). What is the coefficient for age for the 1,000th gene? Make a plot of the data and fitted values for this gene. Does the model fit well?

BiocManager::install("limma")
library(limma)

table(pdata_bm$age)
is.na(pdata_bm$age)

p_data_bm_filtered = subset(pdata_bm, !is.na(age))

edata_no_na = subset(edata,  select=p_data_bm_filtered[,1])

design = model.matrix(~ age, data = p_data_bm_filtered)



print(design)
fit = lmFit(edata_no_na, design)

age_coefficients = fit$coefficients[, "age"]
coefficient_1000 = age_coefficients[1000]
coefficient_1000

# Get the fitted values
fitted_values = fit$fitted.values[1000, ]

# Plotting
plot(age_data$age, edata_sub[1000, ], 
     xlab = "Age", ylab = "Expression of Gene 1000", 
     main = "Gene 1000: Expression vs Age")
points(age_data$age, fitted_values, col = "red", pch = 19)
legend("topright", legend = "Fitted values", col = "red", pch = 19)


## Question 8 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)


p_data_bm = subset(pdata_bm, !is.na(age))

edata = subset(edata,  select=p_data_bm[,1])

mod_adj = model.matrix(~ pdata_bm$age + as.factor(pdata_bm$tissue.type))
fit_limma = lmFit(edata,mod_adj)
fit_limma$coefficients[1,]

### Question 10 

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

set.seed(33353)


library(devtools)
library(Biobase)
#install.packages("sva")
library(sva)
library(bladderbatch)
library(snpStats)

edata = log2(edata +1)
edata = edata(rowMeans(edata > 1))
