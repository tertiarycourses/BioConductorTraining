#biocLite("AnnotationHub", lib=bio)
#biocLite("rtracklayer", lib=bio)


library(AnnotationHub)
library(rtracklayer)

ah=AnnotationHub()

human=subset(ah, species == "Homo sapiens")  # extract subset only human data
human

# function masked by other package, must restart R
query(human, "epidermal")                   # look at EGFR
query(human, c("epidermal","H3K4me3"))      # look at histone

humanEpi=query(ah, c("epidermal","H3K4me3"))
humanEpi

humanEpi[1]
epiBroadPeak=humanEpi[[1]]   #take qst dataset a GRanges (broadpeak)

humanEpi[2]
epiNarrowPeak=humanEpi[[2]]                               #(narrowpeak)

# we compare broadpeak and narrowpeak

epiBroadPeak		#look at the GRanges


summary(width(epiBroadPeak))    
table(width(epiBroadPeak))        

summary(width(epiNarrowPeak))
table(width(epiNarrowPeak))         # data processing has forced most peaks abt 150

peaks=epiNarrowPeak

################ Cross Reference with Genes ################

qhs=query(ah, "RefSeq")   # genes of many species
qhs
qhs$genome

qhuman=subset(qhs, species == "Homo sapiens") #extract only human genes
qhuman

genes=qhs[[1]]   # download first dataset
genes

table(table(genes$name))    # most genes have single transcipt

prom=promoters(genes)
table(width(prom))

args(promoters)    # how we get promoters, arguments for fxn


findOverlaps(prom, peaks)  # how many promoters overlap peaks
ov=length(findOverlaps(prom, peaks))
ov
