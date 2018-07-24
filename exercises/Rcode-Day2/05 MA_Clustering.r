exdat=bodymap.eset   # just click the file
pdata=pData(exdat)
edata=as.data.frame(exprs(exdat))
fdata=fData(exdat)

#################################################################################
#	Cleaning the data
#################################################################################

### check the data for skewed data###
hist(edata[,1],col=2)   # highly skewed data
hist(edata[,1],col=2,breaks=1000) # adjust the number of breaks

### transform the data###
hist(log2(edata[,1]+1),col=2,breaks=100)
# transform the data into log2 (add 1) since we get -Inf since there are zero values in the dataset
hist(log2(edata[,1]+1),col=2,breaks=100,xlim=c(1,15),ylim=c(0,400))
# we choose to zoom in and specify the range we want to see --> remove the zero (high values)

### get rid of the zeros in the data (low values)
hist(rowSums(edata==0),col=2)
# we see that a rows that have values equal to zero
low_genes=rowMeans(edata)<5  # expression less than 5 is considered low
table(low_genes)
filt_edata=filter(as.data.frame(edata),!low_genes) # filter off the low genes
summary(edata)

low_genes2=rowMedians(as.matrix(edata))<5  # use medians instead of means to see if any diff
table(low_genes2)
table(low_genes2,low_genes) # we see most of them have low means and medians

filt_edata2=filter(edata, !low_genes2)


dim(filt_edata2) # compare the two datasets
dim(filt_edata)

#######################################################################################
#	Clustering
########################################################################################

edata=edata[rowMeans(edata)>5000,]
edata=log2(edata+1)


dist1=dist(t(edata)) # transpose the data and calculate distances

 
heatmap(as.matrix(dist1),ColV=NA,RowV=NA)   # visualize the distances via a heatmap

####################################################################
# herarhical clustering
###################################################################
library(dendextend)


hclust1=hclust(dist1)
plot(hclust1)
plot(hclust1, hang=-1) # put all the labels below
dend=as.dendrogram(hclust1)
dend=color_labels(hclust1,4,1:4)
# we put the dendogram as 4 clusters and we colour them differently --> can change "4" to "3" and so on
plot(dend)

labels_colors(dend)=c(rep(1,10),rep(2,9)) # define the range to colour diffrenetly
plot(dend)

#######################################################################
# kmeans clustering
#######################################################################
kmeans1=kmeans(edata, centers=3)
matplot(t(kmeans1$centers),col=1:3,type="l",lwd=3)
# can see one cluster that has high values
# each genes compared to cluster centers and assigned to closest to it
table(kmeans1$cluster)

newdata=as.matrix(edata)[order(kmeans1$cluster),] # putting the data in a matrix first
dim(newdata)
heatmap(newdata,ColV=NA, RowV=NA) # error > file too big
# those which are similiar are clustered together --> can see in heatmap

#clust1hm=as.matrix(edata)[which(kmeans1$cluster==1), ] # too big file
#dim(clust1hm)
#heatmap(clust1hm,ColV=NA, RowV=NA)

clust2hm=as.matrix(edata)[which(kmeans1$cluster==2), ]
dim(clust2hm)
heatmap(clust2hm,ColV=NA, RowV=NA)

clust3hm=as.matrix(edata)[which(kmeans1$cluster==3), ]
dim(clust3hm)
heatmap(clust3hm,ColV=NA, RowV=NA)