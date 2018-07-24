source("http://www.bioconductor.org/biocLite.R")
bio="C://Users//user//Documents//R//win-library//3.4//bioconductor"


#### change environment
R.home() # see where R is installed
.libPaths()

# add to the library paths (temporarily)
.libPaths(c(.libPaths(), bio)) 
.libPaths()

#biocLite()  
biocLite(lib=bio)

############################# Genetic GUIs ######################################

# affylmGUI/oneChannelGUI - GUI for affymetrix data
library(affylmGUI)
affylmGUI()
library(oneChannelGUI)
oneChannelGUI()

# maGUI - GUI for Microarray Data Analysis and Annotation
# library(maGUI)
# maGUI()


# isogeneGUI - A graphical user interface to conduct a dose-response analysis of microarray data
library(IsoGeneGUI)
IsoGeneGUI()


# neaGUI - An R package to perform the network enrichment analysis (NEA).
library(neaGUI)
neaGUI()


# limmaGUI-  facilitate analysis of microarrays and miRNA/RNA-seq data
library(limmaGUI)



# OLINgui - OLIN: optimized normalization, visualization and quality testing of two-channel microarray data.
library(OLINgui)
OLINgui()


# RNASeqGUI-identification of differentially expressed genes
library(RNASeqGUI)
RNASeqGUI()


# SeqGrapheR - visualize cluster of DNA seq read
# library(SeqGrapheR)    > rggobi error
# SeqGrapheR()
