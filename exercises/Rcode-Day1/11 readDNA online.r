library(GenomicRanges)
library(Biostrings)

### downloaded data from database in fasta format
browseURL("https://www.ensembl.org/info/data/ftp/index.html")


### get NCBI files
#biocLite("reutils", lib=bio)

library(reutils)
gene1=efetch("NM_009790", db="nucleotide", rettype="fasta")  # get file from NCBI
gene1

# p53 mouse
gene2=efetch("AB021961", db="nucleotide", rettype="fasta")  # get file from NCBI
gene2

# SNP
gene3=efetch("NC_001477", db="nucleotide", rettype="fasta")  # get file from NCBI
gene3


data=egquery("mouse[orgn]")
data

### get pubmed files
#biocLite("easyPubMed", lib=dir)
library(easyPubMed)
pubmedid=get_pubmed_ids("brca1")
pubmedid

pubmedxml=fetch_pubmed_data(pubmedid)
pubmedxml

pubmedlist=articles_to_list(pubmedxml)
length(pubmedlist)
pubmedlist[1]

pubmedArt=article_to_df(pubmedlist[1])

#### create a a loop

pubmedArtALL=article_to_df(pubmedlist[1])

for(i in 2:length(pubmedlist)){
  pubmedArt=article_to_df(pubmedlist[i])
  pubmedArtALL=rbind.data.frame(pubmedArtALL,pubmedArt)
  }

nrow(pubmedArtALL)
View(pubmedArtALL)

pubmeddwl=batch_pubmed_download("brca1", batch_size=3) #batch_size: maximum number of records to be saved in a single xml files.
pubmeddwl

### get Entrez files
#biocLite("rentrez", lib=bio)
library(rentrez)
gene3=entrez_fetch(db="nucleotide", id=167843256, rettype="fasta")
gene3

snp1=entrez_fetch(db="nucleotide", id=23006943, rettype="fasta")
snp1


### parsing the fasta file
#biocLite("seqinr", lib=bio)
library(seqinr)

# read downloaded files > convert format
snp1.1=read.fasta(file=textConnection(snp1), as.string=T)
snp1.1

# read fasta file from local system
myseq=read.fasta(file="C://Users//user//Desktop//tertiary//dataSets//BIOdatasets//AK292608.fasta")

### write a seq into Fasta
write.fasta(names="SNP-1", sequences=snp1, file.out="snp1.fasta") # check your MyDocuments file

##### DNA analysis

class(snp1)
class(snp1.1)
snp111=strsplit(snp1.1$NZ_AAAP01001546.1,"")
mysequence=snp111[[1]]
mysequence

length(mysequence)
table(mysequence)
GC(mysequence)

count(mysequence,1)    # gives us nucleotide count
count(mysequence,2)    # gives us diclueotide count
count(mysequence,3)    # gives us triclueotide count

translate(mysequence)  # translate into amino acid sequence
tt=paste(translate(mysequence),collapse="")

############################ CHALLENGE #########################

# perform DNA analysis on gene1 and myseq (above)
#FTO gene