library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(GenomicRanges)
library(GenomicAlignments)

#RNA sequencing experiment in the passilaBamSubset package

# our 2 files of interest
fl1=untreated1_chr4()  
fl2=untreated3_chr4()

# coverage function to tally up the base pair coverage
x=readGAlignments(fl1)
xcov=coverage(x)
xcov

# extract one element of the list
xcov$chr4

# zoom in now to range which is near this gene of interest, LGS.
z <- GRanges("chr4",IRanges(456500,466000))
xcov[z]#subset the coverage of the region of interest

par(mfrow=c(1,1))

#Plot of the coverage arround the region of interest
xnum <- as.numeric(xcov$chr4[ranges(z)])#Uncompress the coverage
plot(xnum)


# do the same for another file: 
# we now have pairs of reads which corresponds to one fragment.
y=readGAlignments(fl2)
ycov <- coverage(y)
ynum <- as.numeric(ycov$chr4[ranges(z)])
plot(xnum, type="l", col="blue", lwd=2)
lines(ynum, col="red", lwd=2)

#zoom in on a single exon, between the area of 6000 base pairs:
  
plot(xnum, type="l", col="blue", lwd=2, xlim=c(6200,6600))
lines(ynum, col="red", lwd=2)

################### extracting gene of interest:
# extract lgs gene from the transcript database TxDb.Dmelanogaster.UCSC.dm3.ensGene 
# on Bioconductor, but first we need to look up the Ensembl gene name.

library(biomaRt)
#load the drosophila ensemble gene BioMart.
m <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
lf <- listFilters(m)
lf[grep("name", lf$description, ignore.case=TRUE),]

#get the ensembl gene name
map <- getBM(mart = m,
             attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = "external_gene_name", 
             values = "lgs")
map

#Now we extract the exons for each gene, and then the exons for the gene lgs.
#get the exons out of the transcript database.
library(GenomicFeatures)
grl <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by="gene")
gene <- grl[[map$ensembl_gene_id[1]]]
#View the 6 exons of lgs gene
gene


#Plot each exon as an arrow
rg <- range(gene)

plot(c(start(rg), end(rg)), c(0,0), type="n", xlab=seqnames(gene)[1], ylab="")
arrows(start(gene),rep(0,length(gene)),
       end(gene),rep(0,length(gene)),
       lwd=3, length=.1)




# the lgs gene is on the minus strand. 
# add a line which corrects for minus strand genes:
# If it’s a plus strand, then use code=2 and that means put 
# the arrow head at the end. If it’s a minus strand gene, use code=1 and that means to put an arrow at the start.

plot(c(start(rg), end(rg)), c(0,0), type="n", xlab=seqnames(gene)[1], ylab="")
arrows(start(gene),rep(0,length(gene)),
       end(gene),rep(0,length(gene)),
       lwd=3, length=.1, 
       code=ifelse(as.character(strand(gene)[1]) == "+", 2, 1))


######################## visualizing genomic data: 

#i will show here how to make the coverage plots as before:

library(Gviz)
#set up a Genome Axis Track
gtrack <- GenomeAxisTrack()
#specify an Annotation Track
atrack <- AnnotationTrack(gene, name = "Gene Model")
plotTracks(list(gtrack, atrack))

#Convert coverage to GRanges object
xgr <- as(xcov, "GRanges")
ygr <- as(ycov, "GRanges")
#we have zero coverage for the first four chromosomes
#chromosome four starts out with 891 base pairs of zero coverage.
xgr
ygr

#create two datatracks
#plot coverage which overlap z
dtrack1 <- DataTrack(xgr[xgr %over% z], name = "sample 1")
dtrack2 <- DataTrack(ygr[ygr %over% z], name = "sample 2")
plotTracks(list(gtrack, atrack, dtrack1, dtrack2))
plotTracks(list(gtrack, atrack, dtrack1, dtrack2), type="polygon")

# using ggbio
library(ggbio)
#autoplot gene model
autoplot(gene)
autoplot(fl1, which=z)
autoplot(fl2, which=z)

