##### Using GRanges with genes

library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(ERBS)
library(Homo.sapiens)

#browseVignettes("GenomicRanges")   # read about tools for GenomicRanges

ghs=genes(Homo.sapiens)
ghs

# isolate those on the chrX (chromosome of interest)
ghs[seqnames(ghs)=="chrX"]


ov1=findOverlaps(HumanGenome1, ghs)   #HumanGenome1 from AnnotationHub.r
ov1

#confirm if they overlap (the first hit)
HumanGenome1[1, ]
ghs[20763, ]


ov2=findOverlaps(HumanGenome2, ghs)
ov2


index= queryHits(ov1)     # only want queries that overlap
egfrOV=HumanGenome1[index,] #subset EGFR overlap parts
egfrOV                                    

egfrOVgr= granges(egfrOV)			# we only extract the region data
egfrOVgr


# precede will show genes closest but precedes our site
# precede aslo recognizes if it is "+" or "-" strand
res= precede(egfrOVgr, ghs)
res


ghs[index[1:3]]
egfrOVgr[1:3]

tssgr = resize(ghs, 1)  # reduce all ranges to just 1bp > just want start position
tssgr

d=distanceToNearest(egfrOVgr, tssgr)
queryHits(d)
dist= values(d)$distance

hist(dist,nc=1000, xlim=c(0,100000))  #visual histogram


index=subjectHits(d)[dist<1000]
hist(index)

#### Annotating Genes

keys = as.character(values(tssgr[index])$GENEID)
tssgr[index, ]
keytypes(Homo.sapiens)
columns(Homo.sapiens)
res=select(Homo.sapiens, keys=keys, keytype="GENEID", columns="SYMBOL")
res[1:10, ]


############

df <- data.frame(chr="chr1", start=1, end=1,
                 strand=c("+","-","+","*","."),
                 )

gr= makeGRangesFromDataFrame(df)