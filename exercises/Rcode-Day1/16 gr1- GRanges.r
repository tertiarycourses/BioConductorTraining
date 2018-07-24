# Genomic Ranges package Genomic ranges

#biocLite("GenomicRanges", lib=bio)

library(GenomicRanges)
??GRanges


gr1=GRanges(seqnames=c("chr1"), strand=c("+","-","+"), ranges=IRanges(start=c(1,3,5), width=3))
gr1

gr2=GRanges("chrZ",IRanges(start=c(5,10),end=c(35,45)), strand="+", seqlengths=c(chrZ=100L))
# we specify that chrZ is 100bp long
gr2

seqinfo(gr2)   # info about the genome > 100bp long

### GRangesList   --> groups GRanges together

gr1=GRanges(seqnames = "chr1", ranges = IRanges(start=1:4,width=3))
gr2=GRanges(seqnames = "chr2", ranges = IRanges(start=1:4,width=3))
grL=GRangesList(gr1,gr2)
grL

length(grL)    # length of GRanges
grL[[1]]       # get the first GRange

mcols(grL)$value = c(5,7)
grL

grL[1]      #subset, get first Grange
grL[[1]]    
#grL$gr1     # subset, get first GRange --> works just like list in BaseR


start(grL)    # we get an integer list, we get the start
seqnames(grL) 

elementNROWS(grL)   # how long is the elements in each Grange

shift(grL, 10)        # shift entire GRange by 10

### functions in GRanges

shift(gr2,10)
shift(gr2, 80)  # contain value outside sequence bounds
trim(shift(gr2, 80))   # prevent going outside range

flank(gr1, 5)  # flanking sequence --> relative to direction of transcription
promoters(gr1)  # 200 bases are downstream, 2000 bases upstream

gaps(gr1) # all parts of chromosomes not covered by GRanges and IRanges

#### findOverlaps function

gr1=GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5))
gr2=GRanges("chrZ",IRanges(c(19,33),c(38,35)))
gr1
gr2

fo=findOverlaps(gr1,gr2)
fo

### GRanges metadata and adding more data

mcols(gr2)
mcols(gr2)$value=c(-1,4)   # adding extra column to GRange
# can add more information to our GRanges for more details on DNA

mcols(gr2)$value=NULL    # removes out added columns


seqinfo(gr1)    

seqlengths(gr1)=c("chr1"=10)
seqinfo(gr1)

seqlevels(gr1)

seqlevels(gr1)=c("chr1", "chr2")           # set the levels
seqnames(gr1)=c("chr1", "chr2", "chr1")

seqlevels(gr1)

seqlevels(gr1)=c("chr2", "chr1") # allows sorting of chromosomes
sort(gr1)

genome(gr1)= "hg19"
gr1
seqinfo(gr1)

grX=gr1  # make a copy of gr
genome(grX)="hg18"


##### summary
seqinfo()
seqlevels()
seqnames()
seqlenghts()
genome()


gr11=GRanges(seqnames = "chr11", ranges = IRanges(start=1:4,width=30))
gr22=GRanges(seqnames = "chr22", ranges = IRanges(start=5:8,width=60))
grL=GRangesList(gr11,gr22)
grL

seqinfo(grL)
seqlevels(grL)
seqnames(grL)
seqlengths(grL)=c("chr11"=10,"chr22"=100)
genome(grL)="human"
