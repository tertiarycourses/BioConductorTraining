# rtracklayer is package for importing data

#biocLite("rtracklayer", lib=bio)
library(rtracklayer)

library(AnnotationHub)
ahub=AnnotationHub()
ahub

table(ahub$rdataclass)

### import the BigWigFile data type
ahub.bw = subset(ahub, rdataclass == "BigWigFile" & species=="Homo sapiens") # get only the homosapiens BigWig files
ahub.bw[1]
bw = ahub.bw[[1]]   # extracting the first file

bw

gr.chr22=import(bw, which=GRanges("ch22", ranges = IRanges(1, 10^8)))   # put into memory this range
# we import data and get 1.3 million GRanges

rle.chr22=import(bw, which=GRanges("ch22", ranges = IRanges(1, 10^8)), as = "Rle")
# we import data and get a Rle list

rle.chr22$chr22   # get only the chromosome 22 info


### import the ChainFile data type
ahub.chain = subset(ahub, rdataclass=="ChainFile" & species=="Homo sapiens")
query(ahub.chain, c("hg18","hg19"))

chain=query(ahub.chain, c("hg18", "hg19"))
chain[1]   # get the metadata
chain[[1]] # download the data

#gr.hg18=liftOver(gr.chr22, chain)
#class(gr.hg18)  # converted to Granges

#length(gr.hg18)
#length(gr.chr22)

#table(elementNROWS(gr.hg18))
