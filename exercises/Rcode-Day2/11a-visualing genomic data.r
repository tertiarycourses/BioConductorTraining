
######################### ggbio package #####################

#(devtools)
#install_github("genomicsclass/ERBS", lib=dir)

#biocLite("ERBS", lib=dir)
??ERBS
library(ERBS)
data(HepG2)
HepG2   # this is Grange

#biocLite("GenomeInfoDb", lib=dir)
library(GenomeInfoDb)
seqlevels(HepG2, force=TRUE) = paste0("chr", 1:22)  # isolate only the autosomes
HepG2

#biocLite("ggbio", lib=dir)
library(ggbio)
autoplot(HepG2, layout="karyogram")


### repeat for another data set
data(GM12878)
seqlevels(GM12878, force=TRUE) = paste0("chr", 1:22)
autoplot(GM12878, layout="karyogram")

# karyogram 
p.ideo = Ideogram(genome = "hg19") 
p.ideo  #chr 1 is automatically drawn by default (subchr="chr1") 

library(GenomicRanges) 
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))
#Highlights a region on "chr2"


data(hg19IdeogramCyto, package = "biovizBase")
head(hg19IdeogramCyto)

autoplot(hg19IdeogramCyto, layout = "karyogram", cytoband = TRUE)


### more examples
set.seed(1)
N <- 1000
library(GenomicRanges)
gr <- GRanges(seqnames = 
                sample(c("chr1", "chr2", "chr3"),
                       size = N, replace = TRUE),
              IRanges(
                start = sample(1:300, size = N, replace = TRUE),
                width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N, 
                              replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"), 
                              size = N, replace = TRUE),
              pair = sample(letters, size = N, 
                            replace = TRUE))
idx <- sample(1:length(gr), size = 50)


set.seed(123)
gr.b <- GRanges(seqnames = "chr1", IRanges(start = seq(1, 100, by = 10),
                                           width = sample(4:9, size = 10, replace = TRUE)),
                score = rnorm(10, 10, 3), value = runif(10, 1, 100))
gr.b2 <- GRanges(seqnames = "chr2", IRanges(start = seq(1, 100, by = 10),
                                            width = sample(4:9, size = 10, replace = TRUE)),
                 score = rnorm(10, 10, 3), value = runif(10, 1, 100))
gr.b <- c(gr.b, gr.b2)

###  gr-group

gra <- GRanges("chr1", IRanges(c(1,7,20), end = c(4,9,30)), group = c("a", "a", "b"))
## if you desn't specify group, then group based on stepping levels, and gaps are computed without
## considering extra group method
p1 <- autoplot(gra, aes(fill = group), geom = "alignment")
## when use group method, gaps only computed for grouped intervals.
## default is group.selfish = TRUE, each group keep one row.
## in this way, group labels could be shown as y axis.
p2 <- autoplot(gra, aes(fill = group, group = group), geom = "alignment")

tracks('non-group' = p1,'group.selfish = TRUE' = p2)

### gr-facet strand

autoplot(gr, stat = "coverage", geom = "area", 
         facets = strand ~ seqnames, aes(fill = strand))

### gr - Circle layout


seqlengths(gr) <- c(400, 500, 700)
values(gr)$to.gr <- gr[sample(1:length(gr), size = length(gr))]
idx <- sample(1:length(gr), size = 50)
gr <- gr[idx]
ggplot() + layout_circle(gr, geom = "ideo", fill = "gray70", radius = 7, trackWidth = 3) +
  layout_circle(gr, geom = "bar", radius = 10, trackWidth = 4, 
                aes(fill = score, y = score)) +
  layout_circle(gr, geom = "point", color = "red", radius = 14,
                trackWidth = 3, grid = TRUE, aes(y = score)) +
  layout_circle(gr, geom = "link", linked.to = "to.gr", radius = 6, trackWidth = 1)



######### Interval ranges

set.seed(1)
N <- 100
ir <-  IRanges(start = sample(1:300, size = N, replace = TRUE),
               width = sample(70:75, size = N,replace = TRUE))
## add meta data 
df <- DataFrame(value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
                sample = sample(c("Normal", "Tumor"), 
                                size = N, replace = TRUE),
                pair = sample(letters, size = N, 
                              replace = TRUE))
values(ir) <- df
ir

p1 <- autoplot(ir)
p2 <- autoplot(ir, aes(fill = pair)) + theme(legend.position = "none")
p3 <- autoplot(ir, stat = "coverage", geom = "line", facets = sample ~. )
p4 <- autoplot(ir, stat = "reduce")
tracks(p1, p2, p3, p4)


############ Genomic ranges

set.seed(1)
N <- 100
gr <- GRanges(seqnames = 
                sample(c("chr1", "chr2", "chr3"),
                       size = N, replace = TRUE),
              IRanges(
                start = sample(1:300, size = N, replace = TRUE),
                width = sample(30:40, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N, 
                              replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"), 
                              size = N, replace = TRUE),
              pair = sample(letters, size = N, 
                            replace = TRUE))
grl <- split(gr, values(gr)$pair)


## default gap.geom is 'chevron'
p1 <- autoplot(grl, group.selfish = TRUE)
p2 <- autoplot(grl, group.selfish = TRUE, main.geom = "arrowrect", gap.geom = "segment")
tracks(p1, p2)


################### Gvis package ###################

### genome track plotting

library(Gviz)
data(cpgIslands)
chr="chr7"
genome="hg19"

### we add tracks layer by layer

# simple annotation track
atrack=AnnotationTrack(cpgIslands, name="CpG")   # plotTracks(atrack)
# gtrack represents genomic coordinates
gtrack=GenomeAxisTrack()                        # plotTracks(list(atrack,gtrack))
# ideogram provides overall orientation              
itrack=IdeogramTrack(genome=genome, chromosome=chr) #plotTracks(list(atrack,gtrack,itrack))

data(geneModels)
grtrack=GeneRegionTrack(geneModels, genome=genome, chromosome=chr, name="Gene Model")

tracks=list(itrack, gtrack, atrack, grtrack)
plotTracks(tracks)

### zooming in and out
plotTracks(tracks, from=2.5e7, to=2.8e7)
### adding sequence data
library(BSgenome.Hsapiens.UCSC.hg19)
strack=SequenceTrack(Hsapiens, chromosome=chr)

plotTracks(c(tracks, strack), from=26450430, to=26450490, cex=.8)

#### plotting individual gene ESRRA

#biocLite("Gviz", lib=dir)
#biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene", lib=dir)

library(ERBS)
library(Gviz)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


txdb=TxDb.Hsapiens.UCSC.hg19.knownGene


library(Homo.sapiens)
eid=select(Homo.sapiens, keys="ESRRA", keytype="SYMBOL", columns="ENTREZID")
eid    # we get identity the gene ESRRA 

allg=genes(txdb)
esrraAddr=genes(txdb, filter=list(gene_id=2101))  # get the value of 2101 from eid
esrraNeigh = subsetByOverlaps(allg, esrraAddr+500000)
esrraNeigh$symbol=mapIds(Homo.sapiens, keys=esrraNeigh$gene_id, keytype="ENTREZID", column="SYMBOL")
# we obtain addresses for ESRRA gene body and collect the neighbouring genes, get the symbols

plotTracks(GeneRegionTrack(esrraAddr, showId=TRUE))
plotTracks(GeneRegionTrack(esrraNeigh, showId=TRUE))

idxTrack=IdeogramTrack(genome="hg19", chr="chr11")
# allows us to see the chromosome  # this will take about 3 minutes


geneNeigh=GeneRegionTrack(esrraNeigh, showId=TRUE)
plotTracks(list(idxTrack,geneNeigh)) 
# plot our gene of interest and neighbouring genes on chromosome11


############### GenomeGraphs package #################

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

library(GenomeGraphs)
gene= makeGene(id = "NOX1", type = "hgnc_symbol",biomart = human)
transcript=makeTranscript(id="NOX1",type="hgnc_symbol",biomart=human)
ideogram=makeIdeogram(chromosome ="X")
gdPlot(list(ideogram,gene,transcript))


gene2=makeGene(id = "AGT", type = "hgnc_symbol", biomart = human) 
transcript=makeTranscript(id="AGT",type="hgnc_symbol", biomart=human) 
ideogram=makeIdeogram(chromosome ="1") 
gdPlot(list(ideogram,gene2)) 
gdPlot(list(ideogram,gene2,transcript)) 


minbase=230000000 
maxbase=231000000
genesplus=makeGeneRegion(start = minbase, end = maxbase, strand = "+", chromosome = "1",  biomart = human) 
genesminus=makeGeneRegion(start = minbase, end = maxbase, strand = "-", chromosome = "1",  biomart = human) 
genomeAxis=makeGenomeAxis(add53 = TRUE, add35 = TRUE) 

gdPlot(list(ideogram,genesplus,genomeAxis,genesminus,gene))

