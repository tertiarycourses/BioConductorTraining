# Kyto Encyclopedia of genes and genomes
browseURL("http://www.genome.jp/kegg/")

#biocLite("KEGGREST", lib=bio)

library(KEGGREST)
??KEGGREST

listDatabases()    #get a list of all the databases

org=keggList("organism")  # obtain list of organisms
path=keggList("pathway")   # obtain list of pathways
enz=keggList("enzyme")
dis=keggList("disease")


head(org)   # put into a 2 lists
head(path)
head(enz)

querytables=c(listDatabases(), org[1], org[2])
querytables

head(keggList("hsa"))   # ask for every entry of homo sapiens

keggFind("genes", c("EGFR","hsa"))  # search the EGFR gene

egfrKG=keggGet("hsa:3791")
length(egfrKG)
names(egfrKG[[1]])

egfrKG[[1]]$NAME
egfrKG[[1]]$PATHWAY
egfrKG[[1]]$POSITION
egfrKG[[1]]$AASEQ
egfrKG[[1]]$NTSEQ
egfrKG[[1]]$DBLINKS

efgrAA=keggGet("hsa:3791","aaseq")  # get the aminoacid sequence
efgrNT=keggGet("hsa:3791","ntseq")   # get the nucleotide sequence

efgrAA

efgrNT


library(png)
library(grid)

### for your query > look for pathway > gives the ID for the image


efgrpng=keggGet("hsa04144","image")   #get from $position
grid.raster(efgrpng)

# check class(efgrpng) > should be "Array" contains picture
# and not "NULL"

# lets look at other images
brpng=keggGet("hsa05212", "image")
grid.raster(brpng)  # see image of biological pathway of pancreatic cancer, gene realtionships, TF, processes and so on

brpng2=keggGet("hsa05130", "image")
grid.raster(brpng2)  # see image of E.coli infection pathway


##################### Challenge ############################

# perform a query on BRCA1 gene, get pathways images, NT & AA seq