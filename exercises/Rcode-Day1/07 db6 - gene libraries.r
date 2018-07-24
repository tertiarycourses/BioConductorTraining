# metadata - mapping from entre gene identifiers
#biocLite("org.Hs.eg.db", lib=bio)
library(org.Hs.eg.db)   # homo sapiens library
class(org.Hs.eg.db)

org.Hs.eg.db			

keytypes(org.Hs.eg.db)		# identifier tokens

columns(org.Hs.eg.db)

# seach for pubMed IDs related to EFGR gene

select(org.Hs.eg.db, keys="EGFR", keytype="SYMBOL", columns="PMID")

select(org.Hs.eg.db, keys="EGFR", keytype="SYMBOL", columns="GO")
# MF - metabolic function  BP- biological pathway
# CC- cell cycle

select(org.Hs.eg.db, keys="MAP2K4", keytype="SYMBOL", columns="GENENAME")

# ENZYME NOMENCALTURE EC 1.2.1.3
select(org.Hs.eg.db, keys="1.2.1.3", keytype="ENZYME", columns="GENENAME")

############################### CHALLENGE ###############################

# obtain data about

# BRCA1
# LSM1

# Get info about the gene name and other attribute > look at keytypes and columns