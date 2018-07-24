# Rsamtools package library interface with sequencing library, file for SAM and BAM format

#biocLite("Rsamtools", lib=bio)
library(Rsamtools)

bamPath=system.file("extdata","ex1.bam",package="Rsamtools")
bamPath

bamFile=BamFile(bamPath)
bamFile

seqinfo(bamFile)   #bam file info

aln=scanBam(bamFile)  # read the bam file
length(aln)
class(aln)


aln=aln[[1]]       # look at Bam
names(aln)

aln$qname
aln$qwidth
aln$seq

lapply(aln, function(xx) xx[1])


yieldSize(bamFile)=1   # everytime we call scanBam on the bam file we get just 1 read
open(bamFile)
scanBam(bamFile)[[1]]$seq   # repeat this script for different reads, cycles through list
close(bamFile)
yieldSize(bamFile)=NA

#scan bam params
gr=GRanges(seqnames="seq2",ranges=IRanges(start=c(100,1000), end=c(1500,2000)))
gr
params=ScanBamParam(which=gr,what=scanBamWhat())
scanBamWhat()

aln=scanBam(bamFile, param=params)
names(aln)
head(aln[1]$pos)


quickBamFlagSummary(bamFile)


bamView=BamViews(bamPath)  #read collection of files
bamView
aln=scanBam(bamView)
names(aln)


gr=bamRanges(bamView)
aln=scanBam(bamView)
names(aln)
