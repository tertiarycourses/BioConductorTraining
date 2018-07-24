library(seqinr)
library(bio3d)
library(protr)
library(muscle)

myseq <-"MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMKGDYFRYLSEVASGDNKQTTVSNQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTWTSENQGDEGENLYFQ"

myseq <- s2c(myseq)
table(myseq)

###################################################################################
#				PROTEIN STRUCTURE ANALYSIS
###################################################################################



############# handling pdb files


######## installing

#install.packages("bio3d", dep=T)
#install.packages("msa", dep=T)

library(bio3d)
lbio3d()  # list the functions

example(nma)  # execute examples of a particular function
demo("nma")  # see demo

# Ctrl + C to quit and return to R prompt

# saving work
#read.dcd  > save


#################### PDB analysis ############################

browseURL("http://www.pdb.org")


pdb=read.pdb("C://Users//user//Desktop//tertiary//Rprogramming//BioConductor(basic)//Datasets//pdb1jm7.ent")
pdb  # summary of file


library(protr)
library(bio3d)

pdb <- read.pdb("1BG2")

class(pdb)
attributes(pdb)
head(pdb)
head(pdb$atom[, c("x","y","z")])

# C-alpha coordinates in the protein molecule
head(pdb$atom[pdb$calpha, c("resid", "elety", "x","y","z")])

aa321(pdb$seqres)

write.pdb(pdb, file="myPDBfile.pdb")



######### biod3d

library(bio3d)

pdb1 <- read.pdb("1BG2")
pdb2 <- read.pdb("2VVG")
pdb3 <- read.pdb("1MKJ")


s1 <- aa321(pdb1$seqres)
s2 <- aa321(pdb2$seqres)
s3 <- aa321(pdb3$seqres)


########## protein alignment

#raw <- seqbind(seqbind(s1, s2),s3)

#aln <- seqaln(raw, id=c("1BG2","2VVG","1MKJ"),
#exefile='C://Users//user//Desktop//tertiary//Rprogramming//BioConductor(basic)//muscle3.8.31_i86win32.exe')

#aln2html(aln, append=FALSE, file="Myalign.html")

#Now, open the created HTML file in a browser to visualize the results.


######## computing features of protein sequence

pdb1 <- read.pdb("1BG2")
s1 <- aa321(pdb1$seqres)

s1 <- paste(s1, sep="",collapse="")

#amino acid composition 
extractAAC(s1)
protcheck(s1)

#distribution of the hydrophobic and hydrophilic amino acids along its chain
extractAPAAC(s1, props = c("Hydrophobicity", "Hydrophilicity"),
lambda = 30, w = 0.05, customprops = NULL)

#composition descriptor, transition, descriptor, and dipeptide composition
extractCTDC(s1)
extractCTDD(s1)
extractDC(s1)

##################### PLOTTING

tor <- torsion.pdb(pdb)

plot(tor$phi, tor$psi, main="(A) Ramachandran plot 1BG2")

scatter.psi <- tor$psi
scatter.phi <- tor$phi

library(RColorBrewer) # load RColourBrewer package
k <- 10 # define number of colours
my.cols <- rev(brewer.pal(k, "RdYlBu")) # Brew color pallette

smoothScatter(x=scatter.phi, y=scatter.psi,
colramp=colorRampPalette(my.cols), xlim=c(-180,180), ylim=c(-
180,180),xlab="Phi", ylab="Psi", main="(B) Ramachandran plot
1BG2", pch=19, cex=0.00)


###### 3D visualization

library(Rknots)

myprotein <- loadProtein("1BG2")

plotDiagram(myprotein$A, ends = c(), lwd = 2.5)

ramp <- colorRamp(c('blue', ' white', ' red' ))
pal <- rgb(ramp(seq(0,1,length=100)), max=255)
plotKnot3D(myprotein$A, colors=list(pal), lwd=8, radius=0.4,
showNC=TRUE, text=FALSE)
