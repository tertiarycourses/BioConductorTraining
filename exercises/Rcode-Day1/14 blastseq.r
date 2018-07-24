
#########################################################################
#                         NCBI DNA blast
##############################################################################
browseURL("https://www.ncbi.nlm.nih.gov/BLAST/")

# blastSequences(x, database, hitListSize, filter, expect, program, timeout=40, as="data.frame")

#Arguments
#x            >A sequence as a character vector or an integer corresponding to an entrez gene ID. Submit multiple sequences as a length-1 character vector, x = ">ID-1\nACATGCTA\n>ID-2\nAAACCACTT".
#database     > Which NCBI database to use. If not “blastn”, then set as="XML"
#hitListSize  >Number of hits to keep.
#filter      >Sequence filter; “L” for Low Complexity, “R” for Human Repeats, “m” for Mask lookup
#expect      >The BLAST ‘expect’ value above which matches will be returned.
#program     >Which program do you want to use for blast.
#timeout      >Approximate maximum length of time, in seconds, to wait for a result.
#as          >character(1) indicating whether the result from the NCBI server should be parsed to a list of DNAMultipleAlignment instances, represented as a data.frame, or returned as XML.

#biocLite(annotate, lib=bio)

library(Biostrings)
library(seqinr)
library(annotate)

####read fasta file as string 

myseq=readDNAStringSet("C://Users//user//Desktop//tertiary//dataSets//BIOdatasets//AK292608.fasta")

data=blastSequences(myseq, 
                    hitListSize=20, 
                    timeout=40, 
                    program="blastn", 
                    as="data.frame")
names(data)
length(data)

data[1,]

View(data)

#### data bases
# blastn >nucleotide blast  
# blastp > protein blast   
# blastx > translated nucleotide2protein     
# tblastn   protein2translated nucleotide


######### studying the sequence

seq=getSequence(data$Hsp_qseq[[1]]) # access the first sequence
seq[1:20]

### base count
basecount=table(seq)
myseqname=data$Hit_def[[1]]
dotchart(basecount, xlim=c(0, max(basecount)), pch=19, main=paste("Basecount of", myseqname))

#### codon usage of the sequence
codonusage=uco(seq)
dotchart.uco(codonusage, main=paste("Codon usage in", myseqname))
aacount=table(getTrans(seq))
aacount=aacount[order(aacount)]
aacount
names(aacount)=aaa(names(aacount))
dotchart(aacount, pch=19, xlab="Stop and AA counts")
abline(v=1, lty=2)


######## Challenge ##############

# do a BLAST on FTO gene, get the first sequence and analyze as above

library(Biostrings)
fto=DNAString(yoursequence)   #keying in a DNA string