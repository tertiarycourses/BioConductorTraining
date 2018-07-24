library(Biostrings)

s=readDNAStringSet("C://Users//user//Desktop//tertiary//Rprogramming//BioConductor(basic)//Datasets//AK292608.fasta")
seq=toString(s)  # convert to character
seq=unlist(strsplit(seq,""))


basecount1=count(seq,1)    # gives us nucleotide count
basecount2=count(seq,2)    # gives us diclueotide count
basecount3=count(seq,3)    # gives us triclueotide count

myseqname="AK292608"
dotchart(basecount1, xlim=c(0, max(basecount1)), pch=19, main=paste("Basecount of", myseqname))
dotchart(basecount2, xlim=c(0, max(basecount2)), pch=19, main=paste("Dinucleotide of", myseqname))
dotchart(basecount3, xlim=c(0, max(basecount3)), pch=19, main=paste("Codon count of", myseqname))


#### codon usage of the sequence
codonusage=uco(seq)
dotchart.uco(codonusage, main=paste("Codon usage in", myseqname))
aacount=table(getTrans(seq))
aacount=aacount[order(aacount)]
aacount
names(aacount)=aaa(names(aacount))
dotchart(aacount, pch=19, xlab="Stop and AA counts")
abline(v=1, lty=2)


######################### dna online #####################

library(seqinr)
choosebank()
choosebank(infobank=T)[1:4,]
#choosebank("embl") # how to choose a bank

choosebank("genbank")

?query
completeCDS=query("completeCDS", "sp=Arabidopsis thaliana AND t=cds")
# only complete sequences

# To get the names of the 4 first sequences:
sapply(completeCDS$req[1:4], getName)
# To get the 4 first sequences:
sapply(completeCDS$req[1:4], getSequence, as.string = TRUE)


nseq=completeCDS$nelem
nseq

seq=getSequence(completeCDS$req[[1]]) # access the first sequences
seq[1:20]

### base count
basecount=table(seq)
myseqname=getName(completeCDS$req[[1]])
myseqname
dotchart(basecount, xlim=c(0, max(basecount)),pch=19, main=paste("Basecount of", myseqname))

#### codon usage of the sequence
codonusage=uco(seq)
dotchart.uco(codonusage, main=paste("Codon usage in", myseqname))
aacount=table(getTrans(getSequence(completeCDS$req[1])))
aacount=aacount[order(aacount)]
aacount
names(aacount)=aaa(names(aacount))
dotchart(aacount, pch=19, xlab="Stop and AA counts")
abline(v=1, lty=2)


choosebank("swissprot")
laprae=query("leprae", "AC=Q9CD83")   # we know the ascension number
lepraeseq=getSequence(laprae$req[[1]])
ulcreans=query("ulcerans", "AC=A0PQ23")
ulceransseq=getSequence(ulcreans$req[[1]])

choosebank()
dotPlot(lepraeseq, ulceransseq)

      