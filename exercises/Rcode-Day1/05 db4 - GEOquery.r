#GEOquery package to inteface with NCBI Gene expression Omibus 
# repository for public data microarray data
browseURL("https://www.ncbi.nlm.nih.gov/gds/")


#biocLite("GEOquery", lib=bio)

library(Biobase)
library(GEOquery)
??GEOquery


### log in or signup on NCBI website


# S4 class data type
gds858=getGEO("GDS858")   # you must register and sign in into NCBI to access the data #  
# lung epithelial Calu-3 cells infected with the mucoid alginate-producing FRD1 strain
# can check on NCBI website
class(gds858)
gds858
colnames(Table(gds858))
Table(gds858)[1:10,1:6]

#look at the metadata
Meta(gds858)  # can subset individual components
names(Meta(gds858))
met_gds858=Meta(gds858)
names(met_gds858)  # get the metadata available

# list class data type
gse17948=getGEO("GSE17948")   # microarray data on EGFR
class(gse17948)
gse17948

length(gse17948)
names(gse17948)
eData=gse17948[[1]]
eData
class(eData)

names(pData(eData))
exprs(eData)
head(exprs(eData))
pData(eData)
fData(eData)
annotation(eData)

### save into computer desktop
save(gse17948, file="C://Users//user//Desktop//dwi1test.rda")
rm(gse17948)  # remove from memory
gse17948
load('C://Users//user/Desktop//dwi1test.rda')


gse17948p2 = getGEOSuppFiles("GSE17948")  # download the raw data
gse17948p2


### ArrayExpress  -- lots of microarray data

#biocLite("ArrayExpress", lib=dir)
library("ArrayExpress")

# can check ArrayExpress website
browseURL("https://www.ebi.ac.uk/arrayexpress/")

AEset = ArrayExpress("E-MEXP-1416")


