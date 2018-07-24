#biocLite("IRanges")

library(IRanges, lib=bio)
??IRanges


ir1=IRanges(5,10)   # specify the start and end  --> therfore 6bp long
ir2=IRanges(5,width=6)  # same as ir1

start(ir1)
end(ir1)
width(ir1)

ir3=IRanges(start =c(1,3,5), end=c(3,5,7))  # if we know 2 elements, last element can be inferred
ir3
ir4=IRanges(start=c(1,3,5), width=3)   # we get exactly the same as above
ir4

ir=IRanges(start=c(3,5,17), end=c(10,8,20))  # can specify more that one IRange at a time
ir

length(ir)
start(ir)
end(ir)
width(ir)

width(ir3)=1  # reassign width length to 1 in ir3
ir3           # see new values

################## Challenge ###########################

# create an iRange of 5 elements (random start points) 
# at least 10bp long
# then force all the widths to be only 5bp long




names(ir3)=paste("A",1:3, sep="") # add a names column
ir3
length(ir3) # see how many rows
ir3[1]   # subset ir3 to get first row only
ir3["A1"]  # subset by names
ir5=c(ir3,ir4)  # combine IRanges ir3 and ir4
ir5



### Iranges functions

ir1=IRanges(start=c(3,5,17), end=c(10,8,20))  # we create 3 ranges

shift(ir1, -2)    # shift the range by 2bp

narrow(ir1, start=2) # should start this range at the 2nd bp from start

narrow(ir1, end=3) # should end this range at the 3rd bp form start

flank(ir1, width=3, start=TRUE, both=FALSE)   #get left flanking sequence 3bp from the start

flank(ir1, width=3, start=FALSE, both=FALSE)   #get right flanking sequence 3bp from the start

flank(ir1, width=3, start=FALSE, both=TRUE)   # bidirectional sequencing from the start

range(ir1)                      # from beginning to end
reduce(ir1)		       # those bp which are covered by orginal ranges - no gaps
gaps(ir1)		       # gaps in IRanges
disjoin(ir1)		       # set of ranges which has the same coverage as orginal, non overlap,
# contains union of end points of all original ranges

# resize IRanges #

resize(ir1, width=1, fix="start")  # keep the start values and adjust the end values

resize(ir1, width=1, fix="center") # adjust both start and end values


################## Challenge #######################

# create your own IRange of 5 elements at least 20bp long
# at random starts

# repeat the above Iranges commands
# shift, narrow, flank
# range, reduce, gap , join
# resize


### graphical plot of IRanges function


# plotting wih IRanges (original)

plotir=function(ir,i,c="black"){arrows(start(ir)-.5,i,end(ir)+.5,code=3,angle=90,lwd=3, col=c)}
plot(0,0,xlim=c(-3,15),ylim=c(0,8),type="n",xlab="",ylab="",xaxt="n")
axis(1,-5:15)
abline(v=0.30 + .5, col=rgb(0,0,0,0.5))

ir1
plotir(ir1,1,"red")  #plot original IRange for visualization


plotir(shift(ir1,-2),2)
plotir(narrow(ir1,start=2),3)
plotir(narrow(ir1,end=3),4)
plotir(flank(ir1, width=3, start=TRUE, both=FALSE),5)
plotir(flank(ir1, width=3, start=FALSE, both=FALSE),6)
plotir(flank(ir1, width=3, start=TRUE, both=TRUE),7)



### plotting with IRanges (new)

plotRanges=function(x, xlim=x, main=deparse(substitute(x)), 
                    col="black", sep=0.5, ...) {
                    height=1
                    if(is(xlim, "IRanges"))
                        xlim=c(min(start(xlim)), max(end(xlim)))
                        bins=disjointBins(IRanges(start(x), end(x) +1))
                        plot.new()
                        plot.window(xlim, c(0, max(bins)*(height + sep)))
                        ybottom =bins *(sep + height)-height
                        rect(start(x)-0.5, ybottom, end(x)-0.5, ybottom+height, col=col, ...)
                        title(main)
                        axis(1)
                        }

par(mfrow=c(2,1))     # readjust the size of graphics window > maximise it

ir=IRanges(start=c(1,3,7,9),end=c(4,4,8,10))
plotRanges(ir)
plotRanges(reduce(ir))

plotRanges(ir)
plotRanges(disjoin(ir))

plotRanges(ir)
plotRanges(gaps(ir))

# create two new IRanges set
ir1=IRanges(start=c(1,3,5), width=1)
ir2=IRanges(start=c(4,5,6), width=1)

union(ir1,ir2)
intersect(ir1,ir2)

## find overlaps ##

# create two new IRanges set
ir1=IRanges(start=c(1,4,8), end=c(3,7,10))
ir2=IRanges(start=c(3,4), width=2)
ir1
ir2

ov=findOverlaps(ir1,ir2)
ov
# first row, first element --> range1 overlaps range2
# range2 has multiple overlaps

queryHits(ov)
unique(queryHits(ov))

countOverlaps(ir1,ir2) # range1 has 1 overlap and range2 has 2 overlaps

nearest(ir1,ir2) # which IRanges in ir2 are close to ir1
