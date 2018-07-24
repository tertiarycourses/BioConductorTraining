exdat=bodymap.eset   # just click the file
pdata=pData(exdat)
edata=as.data.frame(exprs(exdat))
fdata=fData(exdat)

dim(edata) # 1198 genes from 19 samples
names(edata)

dim(pdata) #19 samples with 6 characteristics
names(pdata)

### investigate differences in gender
table(pdata$gender)
which(pdata$gender=="M") # find all the rows with males in pData

par(mfrow=c(1,2))
boxplot(edata[which(pdata$gender=="M")])
# plot all the males
boxplot(edata[which(pdata$gender=="F")])
# do the same for females


### investigate differences in age
table(pdata$age)
which(pdata$age<50)
which(pdata$age>=50)

par(mfrow=c(1,2))
boxplot(edata[which(pdata$age<50)])
boxplot(edata[which(pdata$age>=50)])

#### genetic expression

names(fdata)
dim(fdata)
fdata$gene[1:10]  # suppose we only want to focus on first 10 genes


#recovering gene names and chromosome positions from ENSEMBL ID

# select the ENSEMBL database with the human dataset
library(biomaRt)
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")

# here's a vector of ENSEMBL gene id's we're interested in (first 10 genes)
ens_ids=featureNames(exdat)[1:10]

# use getBM() to recover other information about genes:
getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position",
                     "end_position", "band"), filters = "ensembl_gene_id", values = ens_ids, mart = ensembl)

dim(fdata) # 52580 genes
str(fdata)
# suppose we focus on the TNMD gene
boxplot(edata[which(fdata$gene=="ENSG00000000005")])
# suppose we focus on the DPM1 gene
boxplot(edata[which(fdata$gene=="ENSG00000000419")])


# suppose we narrow down to NFYA gene from breast samples
boxplot(edata[which(fdata$gene=="ENSG00000000419"),which(pdata$tissue.type=="breast")])
# in this scenario > we only have one sample >. other cases may have more samples

####################################################################################
#                             Expression analysis
####################################################################################
# further analysis
plot(rowMeans(edata))
plot(colMeans(edata))

#isolate highly expressed genes
highGene=subset(edata, rowMeans(edata)>100)  # we can modify this to any value we wish
View(highGene)


# select the ENSEMBL database with the human dataset
library(biomaRt)
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")

# here's a vector of ENSEMBL gene id's we're interested in all the high Genes
ens_ids=rownames(highGene)
ens_ids=as.factor(ens_ids)

# use getBM() to recover other information about genes:
getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position",
                     "end_position", "band"), filters = "ensembl_gene_id", values = ens_ids, mart = ensembl)


dim(highGene)
highGeneM8=as.matrix(highGene)
heatmap(highGeneM8)

require(gplots)
heatmap.2(highGeneM8,RowV=NA,ColV=NA, dendrogram="none", scale="row", trace="none")

#####################################################################################
#                             statistical tests
#####################################################################################

# 2 sample t test
gene1=unlist(edata[which(fdata$gene=="ENSG00000000005")])
gene2=unlist(edata[which(fdata$gene=="ENSG00000000419")])

var.test(gene1,gene2) # check for equal variances
qqnorm(gene1);qqline(gene1) # check for normality
t.test(gene1,gene2, var.equal=TRUE, paired=TRUE)  # since they are the same genes

# correlation
cor(gene1,gene2)

# we take the first 100 genes from breast,brain,heart,liver samples
geneBreast=edata[fdata$gene[1:100],which(pdata$tissue.type=="breast")] 
geneBrain=edata[fdata$gene[1:100],which(pdata$tissue.type=="brain")]
geneHeart=edata[fdata$gene[1:100],which(pdata$tissue.type=="heart")]
geneLiver=edata[fdata$gene[1:100],which(pdata$tissue.type=="liver")]

geneALL=cbind(geneBreast,geneBrain,geneHeart,geneLiver)

cor(geneALL)
pairs(geneALL[,1:4])
# we analyze that when genes are highly expressed in heart they are 
# also highly expressed in liver

#Anova 
# some arrangement to the data
geneBreast=cbind(geneBreast,rep("breast",100))
geneBrain=cbind(geneBrain,rep("brain",100))
geneHeart=cbind(geneHeart,rep("heart",100))
geneLiver=cbind(geneLiver,rep("liver",100))
geneTest=rbind(geneBreast,geneBrain,geneHeart,geneLiver)

colnames(geneTest)=c("eData","sample")
rownames(geneTest)=NULL
str(geneTest)
geneTest=as.data.frame(geneTest)
geneTest$eData=as.numeric(geneTest$eData)
geneTest$sample=as.character(geneTest$sample)
head(geneTest)

boxplot(geneTest$eData~geneTest$sample)
summary(aov(geneTest$eData~geneTest$sample)) 
# > conclusion is no significant difference

# We can create loops for every 100 genes, and repeat for entire dataset of 1000s genes
# and store the pvalues later for analysis to check for difference

#######################################################################################
#				CHALLENGE
#######################################################################################

############## Challenge 1

#using another dataset in the datasets folder 
#> investigate differences in 2 variables of your choice
#> plot the boxplots and visualize the gene expression
#> recover the gene names using Biomart

############## Challenge 2

#> find the high expressing genes
#> using Biomart to find their names
#> plot a heatmap of the gene expression (high genes)

############## Challenge 3
#choose the first 2 genes and perform a t-test

#>choose first 100 genes in 4 tissue types (if applicable)
#>plot the boxplot and visualize the data
#>perform an ANOVA of tissue types