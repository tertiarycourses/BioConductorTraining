#biomaRt  interface with a biomart biological database online
browseURL("http://www.biomart.org")

#biocLite("biomaRt", lib=bio)  # install biomaRt package
library(biomaRt)   # ask you to choose a database(mart) and database
??biomaRt

head(listMarts())   #list datamarts

mart=useMart("ensembl")  # we choose the ENSEMBL database
mart

head(listDatasets(mart))    # list first 6 datasets

ensembl=useDataset("hsapiens_gene_ensembl", mart)    # choose this dataset

# ensembl website
browseURL("http://useast.ensembl.org")

values=c("202763_at","209310_s_at","277500_at")
# list of affymetrix probe ids, we want to get back gene names associated with it


getBM(attributes=c("ensembl_gene_id", "affy_hg_u133_plus_2"), 
      filters="affy_hg_u133_plus_2", 
      values=values, 
      mart = ensembl)

# filters --> affy_hg_u133_plus_2 is the name of the DNAchip
# values --> data associated with these sets
# getBM --> get bio mart

attribute=listAttributes(ensembl)
head(attribute)
nrow(attribute)
View(attribute)

filter=listFilters(ensembl)
head(filter)
nrow(filter)
View(filter)

getBM(attributes=c("ensembl_gene_id",
                   "ensembl_transcript_id",
                   "ensembl_peptide_id"), 
      filters="affy_hg_u133_plus_2", 
      values=values, 
      mart = ensembl)


getBM(attributes=c("go_id","phenotype_description","gene_biotype","affy_hg_u133_plus_2",
                   "start_position","end_position","chromosome_name"),
      filters="affy_hg_u133_plus_2", 
      values=values, 
      mart = ensembl)



getBM(attributes=c("phenotype_description","gene_biotype",
                   "start_position","end_position","chromosome_name"),
      filters="hgnc_symbol",
      values="FTO",
      mart=ensembl)


attributePages(ensembl) #listPages  -> attributes grouped into pages
attributes=listAttributes(ensembl, page="feature_page")
attributes
head(attributes)
nrow(attributes)


###################### exercises with genomes ##########################

######### lookup databases

# databases are queries with get() or mget() > multiple queries

#biocLite("hgu95av.db",lib=bio)
library("hgu95av2.db")
library(help="hgu95av2.db")
ls("package:hgu95av2.db")

mget(c("738_at", "40840_at", "32972_at"), envir=hgu95av2GENENAME)

go=get("738_at", envir=hgu95av2GO)
names(go)

# get("GO:0009117", envir=GOTERM)


get("NOX1", envir=hgu95av2ALIAS2PROBE)


get("F5",hgu95av2ALIAS2PROBE)
get("AGT", hgu95av2ALIAS2PROBE)
get("35245_at",hgu95av2GO)
get("GO:0007155",hgu95av2GO2PROBE)


# reverse a lookup table with revmap()

get("NOX1", envir=revmap(hgu95av2SYMBOL))

get("X",revmap(hgu95av2CHR))   # all the probes for X chromosome
get("12",revmap(hgu95av2CHR))  # all the probes for chr12


############## biomart

library(biomaRt)
listMarts()

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
horse = useMart("ensembl", dataset = "ecaballus_gene_ensembl")
fish = useMart("ensembl", dataset = "gaculeatus_gene_ensembl")


######## getBM() function  
# getBM returns a list of filters
# for selecting genes or SNPs and attributes to return from the
# database.

affyids=c("202763_at", "209310_s_at", "207500_at")
getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol", "chromosome_name",
       "start_position", "end_position", "band"), 
      filters = "affy_hg_u133_plus_2",
        values = affyids, 
      mart = human)

getBM(attributes=c("chromosome_name","start_position","end_position"),
      filters="mgi_symbol",
      values="Cxcl5",
      mart=mouse)

getBM("GYS1", filters="wikigene_name", 
      attributes="chromosome_name",
      mart=horse)


getBM(attributes="go_id",filters="hgnc_symbol",
      values="F5",
      mart=human)

getBM(attributes="hgnc_symbol", 
      filters="go",
      values="GO:0007155",
      mart=human)

getBM(mart=human, attributes=c("band","hgnc_symbol"),
      filters=c("band_start","band_end","chromosome_name"),
      values=list("p21.33","p21.33",6))


values=c("148350")


getBM(mart=human,
      attributes="ensembl_gene_id",
      values=list(148350,148612,8),
      filters=c("start","end","chromosome_name")
)



getBM(mart=human,
      attributes="chromosome_start",
      values=list(148350,148612,TRUE,8),
      filters=c("start","end","with_validated_snp","chromosome_name")
      )


######### attributes and filters

listAttributes(human)
View(listAtributes(human))

listAttributes(mouse)
listAttributes(horse)
listAttributes(fish)


listFilters(human)
listFilters(mouse)
listFilters(horse)
listFilters(fish)


####### getLDS() function 
# getLDS() combines two data marts, for example to homologous
# genes in other species.

getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"),
       filters = "hgnc_symbol", 
       values = "NOX1", 
       mart = human,
       attributesL = c("chromosome_name","start_position","external_gene_name"),
        martL = mouse)

# take note that the mouse gene name is the same as the human one apart from
# capitalisation.


####### getSequence() function 
# The getSequence function looks up DNA or protein sequences by
# chromosome position or gene identifiers

# HGNC option BRCA1

seq=getSequence(id="A2M", 
                type="hgnc_symbol", 
                mart=human, 
                seqType="transcript_exon_intron")
seq

agt=getSequence(id="AGT",
                type="hgnc_symbol", 
                seqType="peptide", 
                mart=human)
agt

### exercise on p53

getBM(mart=human,
      attributes=c("chromosome_name","start_position","end_position"),
      values="TP53",
      filter="hgnc_symbol"
)

getSequence(id="TP53",
            type="hgnc_symbol", 
            seqType="transcript_exon_intron", 
            mart=human)



################# CHALLENGE #########################

# select 3 cancer genes and perform:
# basic getBM and getSequence
# getLDS (human, mouse, horse, mouse) > choose 2 organisms
# getSequence DNA(exon_inton) and peptide sequence