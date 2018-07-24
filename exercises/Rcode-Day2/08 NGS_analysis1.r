reads=readFastq("C://Users//user//Desktop//tertiary//Rprogramming//BioConductor(basic)//Datasets//ERR127302_1_subset.fastq.gz")
sequences=sread(reads)
sequences[1]

myseq1=sequences[1]

sub=myseq1[grep("99.x", names(myseq1))]
length(sub)

#### pattern matching

mypos=matchPattern("TCTGTGTC", myseq1[[1]], max.mismatch=1)   # find patterns match
countPattern("TCTGTGTC", myseq1[[1]], max.mismatch=1)          # count corresponding matches

tmp=c(DNAStringSet("ATGAAGTC"), DNAStringSet(mypos))
consensusMatrix(tmp)

myvpos=vmatchPattern("ATGAAGTC", myseq1, max.mismatch=1)      # find pattern match
Views(myseq1[[1]], start(myvpos[[1]]), end(myvpos[[1]]))       # retrieve results for single entry
sapply(seq(along=myseq1), function(x) as.character(Views(myseq1[[x]], start(myvpos[[x]]), end(myvpos[[x]]))))   # all matches


##### PWM viewing and searching

pwm=PWM(DNAStringSet(c("GCT", "GGT", "GCA")))

library(seqLogo)
seqLogo(t(t(pwm) * 1/colSums(pwm)))

chr=DNAString("AAAGCTAAAGGTAAAGCAAAA")
matchPWM(pwm, chr, min.score=0.9)


myseq2=DNAString(myseq1[[1]])
matchPWM(pwm, myseq2, min.score=0.9)


####################  Challenge #########################

# there are 20k reads in the sequence file > explore the other reads as above
