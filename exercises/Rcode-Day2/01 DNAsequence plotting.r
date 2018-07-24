########################### visualizing DNA ####################################
library(Biostrings)
myseq=readDNAStringSet("C://Users//user//Desktop//tertiary//dataSets//BIOdatasets//AK292608.fasta")
den1=readDNAStringSet("C://Users//user//Desktop//tertiary//dataSets//BIOdatasets//den1.fasta")

width(myseq)
alphabetFrequency(myseq)

barplot(alphabetFrequency(myseq)[1:4],names.arg=colnames(alphabetFrequency(myseq))[1:4])
# plot the alphabet frequncy

letterFrequency(myseq,"TG")  # > letter must be unique
letterFrequency(myseq,"GC") # find the GC content
letterFrequency(myseq,"G")
letterFrequency(myseq,"GTA")

################### this code generates the nucleotide content graph ###########
#x> we define length of each segment here > can be 20 and so on
# alpha> "GC","GA","TG" or any sequence we wish > single/triple letter can as well

dnaPlot=function(x=10, alpha="GC", myseq=myseq){
mylist=list()
mynum=seq(1,width(myseq),x)
        for(i in 1:length(mynum)){
          mylist=append(mylist,substr(myseq,i,i+x))
        }
library(Biostrings)
mylist2=list()
for(i in 1:length(mylist)){
  mylist2=append(mylist2,letterFrequency(DNAString(unlist(mylist[i])),alpha))
                 }
plot(1:length(mylist2),unlist(mylist2), type="l",xlab="segments of x", 
     ylab=c("frequency of",alpha))
}

dnaPlot(200,"GTC",myseq)
dnaPlot(100,"GC",myseq)

#### create a sliding window plot
# this code examines the GC content of the DNA

slidingwindowplotGConly=function(ncl, myseq){
starts <- seq(1, width(myseq)-ncl, by=ncl)
n <- length(starts) # Find the length of the vector "starts"
chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing 
for (i in 1:n) {
library(Biostrings)
chunk <- substr(myseq,starts[i],(starts[i]+ncl-1))
chunk=DNAString(chunk)
chunkGC <- letterFrequency(chunk, "GC")
print(chunkGC)
chunkGCs[i] <- chunkGC                    }
plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}

# we can adjust the range and the accuracy of the plot
slidingwindowplotGConly(ncl=100, myseq=myseq)
slidingwindowplotGConly(ncl=50, myseq=myseq)
slidingwindowplotGConly(ncl=20, myseq=myseq)


############## sliding window plot for nucleotide content

slidingwindowplotATGC=function(ncl, myseq)
{
library(Biostrings)  
starts <- seq(1, width(myseq)-ncl, by = ncl)
n <- length(starts)
chunkGs <- numeric(n)                       
chunkCs <- numeric(n)
chunkAs <- numeric(n)
chunkTs <- numeric(n)
for (i in 1:n) {
  chunk <- substr(myseq,starts[i],(starts[i]+ncl-1))
  chunk=DNAStringSet(chunk)
  chunkGy <- letterFrequency(chunk,"G")
  chunkAy <- letterFrequency(chunk,"A")
  chunkCy <- letterFrequency(chunk,"C")
  chunkTy <- letterFrequency(chunk,"T")
  chunkG=chunkGy/width(chunk)
  chunkT=chunkTy/width(chunk)
  chunkC=chunkCy/width(chunk)
  chunkA=chunkAy/width(chunk) 
 chunkGs[i] <- chunkG
  chunkAs[i] <- chunkA
  chunkCs[i] <- chunkC
  chunkTs[i] <- chunkT
  }
plot(starts,chunkGs,type="b",col="red",ylab="nuleotide content",ylim=c(0,0.65))
points(starts,chunkAs,type="b",col="green")
points(starts,chunkTs,type="b",col="blue")
points(starts,type="b",chunkCs)
}

# we can adjust the range and the accuracy of the plot
slidingwindowplotATGC(ncl=100, myseq=myseq)
slidingwindowplotATGC(ncl=50, myseq=myseq)
slidingwindowplotATGC(ncl=20, myseq=myseq)


################################## Challenge ########################################

# use the above 3 DNA plotting functions to plot and visualize your fasta file
# make sure it is in DNAStringSet format