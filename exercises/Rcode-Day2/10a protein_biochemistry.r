library(ShortRead)
library(seqinr)

hae=read.fasta(file="C://Users//user//Desktop//tertiary//Rprogramming//BioConductor(basic)//Datasets//EGFR(aa).fasta", seqtype="AA")
length(hae)
a1=hae[[1]]
a1
taborder=a1[order(a1)]
names(taborder)=aaa(names(taborder)) # convert the 1 letter to 3 letter ones
dotchart(table(taborder), pch=19, xlab="amino-acid count")

# compute isoelectric points
computePI(a1)

# compute molecular weight
pmw(a1)

# creating hydropahty scores and plot
data(EXP)
names(EXP$KD)=sapply(words(), function(x) translate(s2c(x)))
kdc=EXP$KD[unique(names(EXP$KD))]
kdc
kdc=kdc[order(names(kdc))]

hydro=function(data, coef) {
      f=function(x){
      freq=table(factor(x, levels=names(coef)))/length(x)
      return(coef %*% freq)}
      res=sapply(data,f)
      names(res)=NULL
      return(res)
}

a=hydro(a1,kdc)
aa=aaa(a1)
dat=data.frame(aa,a)

dat
