### Expression set data container

#biocLite("ALL")
library(ALL)
??ALL
?ExpressionSet

data(ALL)
annotation(ALL)
class(ALL)
slotNames(ALL)
experimentData(ALL)
phenoData(ALL)

ALL


exprs(ALL)[1:4, 1:4]

head(sampleNames(ALL))

head(featureNames(ALL))

head(pData(ALL))  #phenotype data info about sample

head(pData(ALL)$sex)
head(ALL$sex)


ALL[,1:5]   # first 5 samples
ALL[1:10,]  # first 10 features
ALL[1:10,1:5]

ids=featureNames(ALL)[1:5]  # id of genes involved in study
ids


####### example 2

#biocLite("hgu95av2.db", lib=bio)   # we get this from annotation(ALL)
library(hgu95av2.db)  # microarray data package
??hgu95av2.db


as.list(hgu95av2ENTREZID[ids])


pD=phenoData(ALL)  #not the same as pData

varLabels(pD)


varLabels(pD)[2]

colnames(pD)[1:3]


#### Summarized Experiment data container

#biocLite("SummarizedExperiment", lib=bio)
#biocLite("airway", lib=bio)
library(airway)
??airway
?SummarizedExperiment


data(airway)
slotNames(airway)
colData(airway)

data(airway)  # sample data set
airway


airway$cell

metadata(airway)

colnames(airway)
head(rownames(airway))

assayNames(airway)
assay(airway, "counts") [1:4, 1:4]  # first 4 genes, first 4 colomns

length(rowRanges(airway))

rowRanges(airway) # GRanges list; each row is a gene, each GRange is the exon of the gene
# almost 64000 genes

sum(elementNROWS(rowRanges(airway))) #almost 750000 exons

start(airway)

gr=GRanges("1",range=IRanges(start=1, end=10^7))
gr

subsetByOverlaps(airway, gr)