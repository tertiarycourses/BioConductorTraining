biocLite("AnnotationHub", lib=bio)

################################## Annotaion Hub database ###############################

library(AnnotationHub)
?AnnotationHub


ah=AnnotationHub()
ah                      # different data resources

ah[1]                   # first dataset metadata

unique(ah$dataprovider)    # resource providers in databases

unique(ah$species)        # look at different species in databases

human=subset(ah, species == "Homo sapiens")  # extract subset only human data
human

query(human, "epidermal")                   # look at EGFR
query(human, c("epidermal","H3K4me3"))      # look at histone

query(human, c("epidermal","H3K4me3"))


humanEpi=query(ah, c("epidermal","H3K4me3"))
humanEpi

human2=display(humanEpi)                        # new display window will appear at sidebar
# select any number of rows of interest and click button "Return rows to R session" #
# this will store the rows in variable human2#

human2

# download that dataset

library('rtracklayer')
dataset=human2[[1]]   # download the dataset

dataset

################################## RefSeq ###########################################

# transcript database (gene annotation)  validated genes

qhs=query(ah, "RefSeq")   # genes of many species
qhs
qhs$genome

qhuman=subset(qhs, species == "Homo sapiens") #extract only human genes
qhuman

qhuman[1]    #hg19 (latest human genome database)
qhuman[2]

HumanGenome1=qhuman[[1]]
HumanGenome1

HumanGenome2=qhuman[[2]]
HumanGenome2


########################## CHALLENGE ##############################################

# Using AnnotationHub >obtain data on Cpg islands on the human genome
# how many islands exist on chromosomes ?
# how many islands exist on chromosome 5 ?

# using hg19 > what is the length of chr 16?

mycpg=query(human, c("CPG","hg19","evoCPG"))
mycpg2=mycpg[1]
mycpg22=mycpg2[[1]]


save(mycpg22, file="C://Users//user//Desktop//mycpg22.rda")

### SNPs

snp=query(human,"snp")
snp2=display(snp)

snp22=snp[[1]]
