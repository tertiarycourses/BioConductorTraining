####################################################################################################
#                          manipulate bioloigcal strings (DNA, RNA, AA)
#####################################################################################################
#biocLite("Biostrings")
library(Biostrings)


####read fasta file as string 
s=readDNAStringSet("C://Users//user//Desktop//tertiary//dataSets//BIOdatasets//AK292608.fasta")

dna1=DNAString("AAGACCTCAGCT")   #keying in a DNA string
dna1

dna2=DNAStringSet(c("ACG","ACT","ATT", "GCA", "GCC")) #creating a set of DNA strings
dna2

# convert from downloaded fasta file from NCBI
snp1=entrez_fetch(db="nucleotide", id=23006943, rettype="fasta")
snp1   # XML format
snp1.1=read.fasta(file=textConnection(snp1), as.string=T)
snp1.1

dna11=DNAString(as.character(snp1.1$NZ_AAAP01001546.1))

# convert from downloaded BLAST RESULT from NCBI
# seq was the result of our blast the first sequence
seq2=paste(seq, collapse="")
dna12=DNAString(seq2)


# randomly generate DNA string
dna3=sample(LETTERS[c(1,3,7,20)], size=1000, replace=TRUE)
dna3 = paste(dna3, collapse="")
dna3=DNAString(dna3)
dna3

rna1=RNAString("GCAUAUUAC")  # RNAString object

prt1=AAString("HCWYHH")      # protein string object

######### inducing mutations

# substitution
dna1[5]="A"
dna1[8]="T"

# addition
dna1=append(dna1, DNAString("T"))
dna1=append(dna1, DNAString(paste(rep("A",20), collapse="")))

dna1[2:4]  # subsequence of the sequence
subseq(dna1, 2,4)  # another method

dna2[1:2]  #subset a DNA string set

dna2[[1]]   # first element of DNA string set

names(dna2)=paste("seq",1:5)  # give names to DNA string set

width(dna2)  #how many bases in string

sort(dna2)   # sort lowest to highest

rev(dna2)   # just reverse the order of the set

rev(dna1)   # reversing the sequence of bases on a string

reverse(dna2)   # true reverse; even in DNA string set

reverseComplement(dna2)   # gives the complement sequence of the DNA set

alphabetFrequency(dna2)   # gives the alpabet frequency of the DNA string set elements

IUPAC_CODE_MAP    # coding map

letterFrequency(dna2, letters = "AC")  # counting occurance of both A and C together in DNA string set
letterFrequency(dna2, letters = "G")
letterFrequency(dna2, letters = "GC")


dinucleotideFrequency(dna2)   # counts the occurance of dinucleotide set in the DNA string
trinucleotideFrequency(dna2)

consensusMatrix(dna2)   # how many strings have particular nucleotides

##########################################################################################################
#				sequence manipulation and alignment
###########################################################################################################
# grep to search matching sequences
sequences1=c(dna1,dna3)
nchar(sequences1) # get length of string

pos1=regexpr("AC", sequences1)
as.numeric(pos1)               # return the match position
attributes(pos1)$match.length  # return the match length

as.numeric(pos1[[1]]) #> to get piece by piece

# change sequence by pattern
gsub("TC", "tc", sequences1)

# extract sequence regions
substring(sequences1,3,10)

# get subsequences
sub1=subseq(dna3,5,8)
sub2=subseq(dna3,5,width=4)
sub3=subseq(dna3,5,-5) # indexing from end

sequences2=DNAStringSet(c("AGGCTCT","AGGTCT","AGGGTGT","AGGCTCT","TTCGGTA"))
duplicated(sequences2) # check for duplicated
unique(sequences2)
sort(sequences2)
order(sequences2)
sequences3=sequences2[2:4]
sequences3=append(sequences2, DNAStringSet("TTCGGTC"))

#############################################################################################################
#                                      pattern matching
#############################################################################################################

library(BSgenome)

available.genomes() # lists all the downloadable genomes in BioConductor website

#biocLite("BSgenome.Scerevisiae.UCSC.sacCer2")
library("BSgenome.Scerevisiae.UCSC.sacCer2")  # getting the yeast genome --> we get an object type
Scerevisiae               # object loaded in library  -> view the sequence data

# nothing is loaded in library

seqnames(Scerevisiae)     # get names of sequence
seqinfo(Scerevisiae)
seqlengths(Scerevisiae)   # get lengths of sequence

Scerevisiae$chr1          # loaded in specific chrmosome sequence into memory

letterFrequency(Scerevisiae$chr1, "GC")   # gives GC content

letterFrequency(Scerevisiae$chr1, "GC", as.prob = TRUE)  # shows GC content as probability

# to find GC content of all the chromosomes
# using bsapply function
# similiar to apply function of base R

param = new("BSParams", X = Scerevisiae, FUN = letterFrequency)  #specify parameter for bsapply
bsapply(param, "GC")   # give argument we wish to put in function
unlist(bsapply(param, "GC"))
unlist(bsapply(param, "GC", as.prob = TRUE)) # express as probability

# to find GC content of yeast genome
sum(unlist(bsapply(param, "GC")))/sum(seqlengths(Scerevisiae))


# create a DNA string
dnaseq = DNAString("ACGTAGCT")
dnaseq

matchPattern(dnaseq, Scerevisiae$chr1)   # finding our sequence on the yeast genome chrI
# retrun object is a view

countPattern(dnaseq, Scerevisiae$chr1)  # count how many matches we have

vmatchPattern(dnaseq, Scerevisiae)   # match against many sets of sequences (entire genome etc)
# we get GRanges

# in this special case for every hit on the forward strand we get a hit on the reverse strand
# this is because our sequence is reverse complement

dnaseq == reverseComplement(dnaseq)   # check if reverse complement


#########################################################################################################################
#						CHALLENGE
#########################################################################################################################

########### Challenge 1

#Obtain data from H3K4me3 histone modification from the H1 cell line from Epigenomic roadmap using AnnotationHub

#Subset these regions to keep only regions mapped to autosomes(chr1 to chr22)

#how many bases do the regions cover?

#repeat this exercise for H3K27me3

########### Challenge 2

#What is the GC content of chr 22 in hg19 (human genome 19)

#CpG islands are clusters of CpGs. What is the observed number of CG dinucleotides for CpG islands in chr 22

#A TATA boc is a DNA element in the form of "TATAAA". Around 25% of genes should have a TATA box in their promoters.
#How many TATA boxes are there on chr22 on hg19 (human genome 19)

############ Challenge 3

# Which of the following sequences are most common on chr11
#   1) ATG   2) TGA   3)TAA  4)TAG