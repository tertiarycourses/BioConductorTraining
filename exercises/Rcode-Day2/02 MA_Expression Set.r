library(gplots)
library(devtools)
library(Biobase)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(dendextend)
library(ALL)

?? ALL

data(ALL)       # summarizing the ALL dataset
annotation(ALL)
class(ALL)
slotNames(ALL)
experimentData(ALL)
phenoData(ALL)

ALL

edata=exprs(ALL)  # extract expression data 
str(edata)
head(edata)
summary(edata)


pdata=pData(ALL)  # extract phenotype data --> describes the samples
str(pdata)
head(pdata)
summary(pdata)
names(pdata)
table(pdata$sex)
table(pdata$sex, pdata$age)


fdata=fData(ALL)  # extract feature data--> describe the genes
str(fdata)
head(fdata)
summary(fdata)
featureNames(ALL)[1:10]

##############################################################################
#  Expression Sets
##############################################################################

load("C://Users//user//Desktop//tertiary//dataSets//BIOdatasets//ExpressionSets//bodymap_eset(human).RData")

bm=bodymap.eset
bm

pdata=pData(bm) # extract phenotype data > describes the samples
dim(pdata)
head(pdata)
str(pdata)        #data exploration
summary(pdata)
table(pdata$gender)
table(pdata$gender, pdata$race)
table(pdata$age)
table(pdata$age, useNA="ifany")  # check if any NAs
sum(pdata==" ")    #another check if any values that are <space>


edata=exprs(bm) #extract expression set
dim(edata)
head(edata)
str(edata)
summary(edata)   # can see highly skewed data

is.na(edata)[1,]   
sum(is.na(edata))   #another check to see if any NAs
gene_na=rowSums(is.na(edata))
table(gene_na)
sample_na=colSums(is.na(edata))
table(sample_na)


fdata=fData(bm) # extract feature data > describes the genes
dim(fdata)
head(fdata)


# to access phenotype table:
phenoData(bm)  #gives information about the table
phenoData(bm)@data  #outputs the table
phenoData(bm)$gender  #gives the gender variable as a vector

# to access gene names:
featureNames(bm)[1:10]  # gives first 10 genes in the count table

# to access count table:
bodymap.count.table <- exprs(bm)
dim(bodymap.count.table)
names(bodymap.count.table)  
as.character(names(bodymap.count.table)) == phenoData(bm)$sample.id  #TRUE (so these are the same thing)
bodymap.count.table[1:10,1:5]  #1st 10 rows and 1st 5 columns of the count table.
