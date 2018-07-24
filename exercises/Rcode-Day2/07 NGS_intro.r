#biocLite('shortRead', lib=bio)

library(ShortRead)

reads=readFastq("C://Users//user//Desktop//tertiary//Rprogramming//BioConductor(basic)//Datasets//ERR127302_1_subset.fastq.gz")
reads

sequences=sread(reads)


################ working with shortread data
sequences[1]
sequences[10000]

substr(sequences,1,5)  # will return a vector for the first 5 bases of each read
table(substr(sequences,1,5))

### filtering data

# suppose we are only interested in sequences beginning with TAGCT (eg PCR primer)

myreads=reads[substr(sequences,1,5)=="TAGCT"]  # sequences,start,stop
myreads  # only 3 results came back

sread(myreads[1])
sread(myreads[2])
sread(myreads[3])


myreadsXX=reads[substr(sequences,65,72)=="TATGCCC"]
myreadsXX


################## filter function
     myfilter=srFilter(function(x){
               substr(sread(x),1,5)=="TAGCT"
               }, name="My Filter")
  
# x is the ShortReadQ object

myreads=reads[myfilter(reads)]


################# matching & counting

dict1=DNAStringSet(sequences)

# suppose we are interested in occurance of a "special sequence" in the reads between pos 22 and 44
dict2=DNAStringSet(substr(sequences,22,44))

hits=vcountPattern("TTGTAATTAAGG", dict1, max.mismatch=1, with.indels=T)
sum(hits)   # number of hits in each read

sread(reads[hits])
id(reads[hits])    # get the ID of our results

############# Challenge #########################

# explore the dataset given -> ERR127302_2_subset.fastq.gz  as above



###################### working shortread package ###################

library(ShortRead)

######## retrive info from fastq files

fastq=readFastq(fastqFile)   # readFasta, writeFasta, writeFastq

seqId=id(fastq)
seqs=sread(fastq)
qualSeq=quality(fastq)
totalReads=length(fastq)

# filter the fastq files

filter=nFilter(threshold=3)   # only reads fewer than 3 Ns

filter2=polynFilter(threshold=20, nuc=c("A","C","T","G"))  # with 20 or more composing of same letters

filter=compose(filter1,filter2)   # combine filters into one

filteredReads=fastq[filter(seqs)]   #apply filter to sequences, remove the bad reads

writeFastq(filteredReads, outputfile)

######## retrive info from bam files

bam=scanBam(bamLoc)[[1]]
names(bam)
scanBamHeader(bamLoc)


cseq=as.character(bam$seq)
cig=bam$cigar
head(cig,2)

qual=bam$qual
head(qual,2)

qname=bam$qname
head(qnae,2)

rname=as.character(bma$rname)
head(rname,2)


## BAM QC
readAligned(bamLoc, type="BAM")

length(aln) # total number of hits

length(unique(id(aln))) #count number of query sequences

length(unique(id([srduplicated(id(aln))])))  # count reads producing multiple hits

table(table(as.character(id(aln)))==1)["TRUE"]  # count reads with unique(single hit) alignment

length(unique(position(aln)))  # count unique positions hit

length(position(aln)[duplicated(position(aln))]))  # count positions with multiple aligning reads


