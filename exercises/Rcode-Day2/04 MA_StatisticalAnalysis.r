load("C://Users//user//Desktop//tertiary//Rprogramming//BioConductor(basic)//Datasets//ExpressionSets//bodymap_eset(human).RData")

exdat=bodymap.eset   # just click the file
pdata=pData(exdat)
edata=as.data.frame(exprs(exdat))
fdata=fData(exdat)

#################### Quality control #####################

### boxplot
boxplot(edata)
boxplot(log(edata))   # log transform
boxplot(log(edata + 1))
        
str(edata[1:100,1:19])

hist(edata[,1]) # all genes for first sample
hist(log(edata[,1]))
hist(log(edata[,1]+1))


### hierarchical clustering

tdata=t(edata)
d=dist(tdata)
plot(hclust(d))

# groups: pdata$gender ; pdata$age; pdata$race

plot(hclust(d), labels=pdata$gender, hang=-1)
plot(hclust(d), labels=pdata$race, hang=-1)
plot(hclust(d), labels=pdata$age, hang=-1)

#### multidimensional scaling

library(MASS)
tdata=t(edata)
d=dist(tdata)
mds=isoMDS(d)
plot(mds$points[,1],mds$points[,2], col=pdata$gender)
plot(mds$points[,1],mds$points[,2], col=pdata$age)
plot(mds$points[,1],mds$points[,2], col=pdata$race)

########################## filtering ###############################
# remove data of poor quality or too low intensity
# non interesting > features not showing any changes

#### fold changes between groups 

group1=edata[which(pdata$gender=="M")]
group2=edata[which(pdata$gender=="F")]

g1=rowMeans(group1)
g2=rowMeans(group2)

par(mfrow=c(2,1))
plot(g1)
plot(g2)

# FoldChange can be calculated by
fc=g1/g2
ind=fc>2   # can change to any value you wish
edata.fc=edata[ind, ]

# extract all data with 2 or more fold change
## draw a herarchical cluster of genes that show at least 2 fold changes between groups

#### filtering by SD
# removing non interesing data

sds=apply(edata,1,sd)
cutsd=quantile(sds, 0.95)
ind=sds>cutsd  # rows that meet our condition > top 5% data
edata.sd=edata[ind,]


####################### clustering ##############################

km=kmeans(edata, centers=5, nstart=10)

plot.km=function(x, cl){
  max.dat=max(x)
  min.dat=min(x)
  par(mfrow=c(ceiling(sqrt(length(cl$size))), ceiling(sqrt(length(cl$size)))))
  for(i in 1:length(cl$size)){
    matplot(t(x[cl$cluster==i,]), type="l", main=paste("cluster:",i), ylab="log expression", col=1, lty=1,
            ylim=c(min.dat, max.dat), pch=".")
  }
}

plot.km(edata,km)

## extract out cluster 5

summary(km)
km

ind=(km$cluster==1)
edata.km=edata[ind, ]


##################### clean and normalize the data #####################

### check the data for skewed data###
hist(edata[,1],col=2)   # highly skewed data
hist(edata[,1],col=2,breaks=1000) # adjust the number of breaks

### transform the data###
hist(log2(edata[,1]+1),col=2,breaks=100)
# transform the data into log2 (add 1) since we get -Inf since there are zero values in the dataset
hist(log2(edata[,1]+1),col=2,breaks=100,xlim=c(1,15),ylim=c(0,400))
# we choose to zoom in and specify the range we want to see --> remove the zero (high values)

### get rid of the zeros in the data (low values)
hist(rowSums(edata==0),col=2)
# we see that a rows that have values equal to zero

plot(rowMeans(edata))
mean(rowMeans(edata))
median(rowMeans(edata))


low_genes=rowMeans(edata)<5  # expression less than 5 is considered low
table(low_genes)

high_genes=rowMeans(edata)>50

filt_edata1=edata[high_genes,]   # filter off the low genes
dim(filt_edata1)  # check the dimensions if R is able to plot

  
###################### creating a model #############################

#biocLite('limma', lib=bio)

library(limma)
design=model.matrix(~pdata$race)  # error in pdata$gender
fit=lmFit(filt_edata1, design)
fit=eBayes(fit)
tp=toptable(fit,2,nrow(filt_edata1))
tp


# first column > row number
# logFC > log transformed fold change between groups
# pvalue > Raw pValue for the feature
# adjPvalue > false discovery rate

### extract significant features

psig=tp$P.Value<0.05
ind=as.numeric(rownames(tp))
ind2=ind[psig]
sign=edata[ind2, ]


####################### plots #########################################

library(limma)

# MA plot
# The plot visualises the differences between measurements taken in two samples, 
# by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values. 
# The MA-plot provides a global view of the differential genes, with the log2 fold change on the y-axis 
# over the mean of normalized counts:
# The MA plot highlights an important property of RNA-seq data. For weakly expressed genes, we have no 
# chance of seeing differential expression, because the low read counts suffer from such high Poisson noise 
# that any biological effect is drowned in the uncertainties from the sampling at a low rate. We can also 
# show this by examining the ratio of small p values (say, less than 0.05) for genes binned by mean normalized 
# count. We will use the results table subjected to the threshold to show what this looks like in a case when 
# there are few tests with small p value.
plotMA(filt_edata1)


# volcano plot
#a type of scatter-plot that is used to quickly identify changes in large data sets composed of 
# replicate data.[1] It plots significance versus fold-change on the y and x axes, respectively.
# A volcano plot combines a measure of statistical significance from a statistical test 
# (e.g., a p value from an ANOVA model) with the magnitude of the change, enabling quick visual 
# identification of those data-points (genes, etc.) that display large magnitude changes that 
# are also statistically significant.
# where the x axis is related to a measure of the strength of a statistical signal, and y axis 
# is related to a measure of the statistical significance of the signal
plot(y=-log10(tp$P.Value), x=tp$logFC, pch=19)
abline(h=2, lty=2)
abline(v=-1, lty=2)
abline(v=1, lty=2)


# heat map
# a graphical representation of data where the individual values contained in a matrix are 
# represented as colors

edata.mat=as.matrix(filt_edata1)
heatmap(edata.mat)
