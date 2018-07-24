# BSgenome package for representing full genomes in BioConductor
#biocLite("BSgenome", lib=bio)

library(BSgenome)
??BSgenome

available.genomes() # lists all the downloadable genomes in BioConductor website
installed.genomes()  # list of available genomes in your computer

#biocLite("BSgenome.Scerevisiae.UCSC.sacCer1", lib=bio) # installing the genomes
library("BSgenome.Scerevisiae.UCSC.sacCer1")  # getting the yeast genome --> we get an object type Scerevisiae               
# object loaded in library  -> view the sequence data

# nothing is loaded in library

seqnames(Scerevisiae)     # get names of sequence

seqlengths(Scerevisiae)   # get lengths of sequence

Scerevisiae$chr1          # loaded in specific chrmosome sequence into memory
class(Scerevisiae$chr1)   # biostrings class

y1=Scerevisiae$chr1
length(y1)

###### human
library("BSgenome.Hsapiens.UCSC.hg19")

seqnames(Hsapiens)

seqlengths(Hsapiens)

hs17=Hsapiens$chr17
length(hs17)