#biocLite("IRanges", lib=bio)

library(IRanges)

ir=IRanges(start=c(3,8,14,15,19,34,40), end=c(14,13,19,29,24,35,46))

library(GenomicRanges)
gir=GRanges(seqnames="chr1", ir)
gir

strand(gir)=c(rep("+",4), rep("-",3))
genome(gir)="hg19"
seqinfo(gir)

#### plotting with GRanges

plotGRanges=function(x, xlim=x, main=deparse(substitute(x)), 
                    col="black", sep=0.5, ...) {
                    height=1
                    if(is(xlim, "GRanges"))
                        xlim=c(min(start(xlim)), max(end(xlim)))
                        bins=disjointBins(IRanges(start(x), end(x) +1))
                        plot.new()
                        plot.window(xlim, c(0, max(bins)*(height + sep)))
                        ybottom =bins *(sep + height)-height
                        rect(start(x)-0.5, ybottom, end(x)-0.5, ybottom+height, col=col, ...)
                        title(main)
                        axis(1)
                        }

plotGRanges(gir, xlim=c(0,60)) # original GRange

# readjust the size of graphics window > maximise it


par(mfrow=c(2,1))

plotGRanges(reduce(gir))  # project all occupied bases into continuous interval
			 # used for reducing complexity
plotGRanges(gir, xlim=c(0,60))


plotGRanges(disjoin(gir))  # set of interval by disjoin of ranges
			  # more complex set of intervals
plotGRanges(gir, xlim=c(0,60))


plotGRanges(gaps(ir), xlim=c(0,60))   # can be area never expressed - exons
plotGRanges(gir, xlim=c(0,60))


#pick out transciption start sites
plotGRanges(resize(gir,1), xlim=c(0,60), col="green") # resize start position to just 1bp upstream
# the last 3 strand are -ve, so start at the end to front
plotGRanges(gir, xlim=c(0,60))

#pick out promoters
plotGRanges(flank(gir,3), xlim=c(0,60), col="purple") # assume each promoter is 3bp upstream
plotGRanges(gir, xlim=c(0,60))

#pick out terminators
plotGRanges(flank(gir,2, start=FALSE), xlim=c(0,60), col="brown")
# downstream end of intervals, 2bp down
plotGRanges(gir, xlim=c(0,60))

