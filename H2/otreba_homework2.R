install.packages('tidyverse')
#install.packages('sva')
#install.packages('patchwork')
#install.packages('RColorBrewer')
install.packages('gplots')
install.packages('devtools')
install.packages('broom')
#install.packages('Biobase')
#install.packages('limma')
#install.packages('edge')
#install.packages('genefilter')
#install.packages('GEOquery')
#install.packages('qvalue')
install.packages('jackstraw')
install.packages('corpcor')

library(corpcor)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase","limma","genefilter","edge","qvalue"))
install.packages('devtools')

library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
library(qvalue)
library(tidyverse)
library(data.table)

load(file='bottomly.Rdata')
ls()

edata <- as.matrix(exprs(bottomly.eset))
dim(edata)

edata <- edata[rowMeans(edata) > 10, ]
edata <- log2(as.matrix(edata) + 1)

library(RColorBrewer)
library(gplots)

my_palette <- colorRampPalette(c("blue", "white", "yellow"))(n = 299)

heatmap.2(edata,
          main = "Heatmap", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "row",
          Colv=FALSE)

dev.off()

#HHHHHHHHHHHHHHHHHHHHHHHHHH
heatmap.2(edata,
          main = "Bottomly et al. Clustered", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="column",     # only draw a row dendrogram
          scale = "column")

#homework no 1!!!!!!!!
#new heatmap, use heatmap.2
#same function with different arguments
#HHHHHHHHHHHHHHHHHHHHHHHHHH

edata<-t(scale(t(edata), scale=FALSE,center=TRUE))
svd.out<-svd(edata)
names(svd.out)

print(paste("Dimension of left singular vectors:",dim(svd.out$u)))
print(paste("Length of singular values:",length(svd.out$d)))
print(paste("Dimension of right singular values:",length(svd.out$v)))

par(mfrow=c(1,1))
plot(svd.out$d,pch=20,ylab="Singular values")
plot(svd.out$d^2/sum(svd.out$d^2)*100,pch=20,ylab="% variance dimension")

plot(1:ncol(edata),svd.out$v[,1],pch=20)

PC<-data.table(svd.out$v,pData(bottomly.eset))

ggplot(PC)+geom_point(aes(x=V1,y=V2,col=as.factor(strain)))

ggplot(PC)+geom_point(aes(x=V1,y=V2,col=as.factor(lane.number)))

ggplot(PC)+geom_point(aes(x=V1,y=V2,col=as.factor(experiment.number)))

#HHHHHHHHHHHHHHHHHHHHHHHHH

#V2 i V3
ggplot(PC)+geom_point(aes(x=V2,y=V3,col=as.factor(strain)))

#homework no 2 !!!!!!
#opisane gdzieś tam

#HHHHHHHHHHHHHHHHHHHHHHHHH





#HHHHHHHHHHHHHHHHHHHHHHHHH

plot(1:ncol(edata),svd.out$u[1,],pch=20, col='red')
#legend(1,1,legend=c('Column 1', 'Column 2'), col=c('red', 'blue'))
points(1:ncol(edata),svd.out$u[2,],pch=20, col='blue')
legend(0.9,0.9,c('Column 1', 'Column 2'), col=c('red', 'blue'))
#plot(1:ncol(edata),svd.out$u[2,],pch=20, col='blue')

#homework no 3 !!!!
#scatter plot for left singular vector

#HHHHHHHHHHHHHHHHHHHHHHHHH

#HHHHHHHHHHHHHHHHHHHHHHHHH

#ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(aes(x=as.factor(strain), y=V1))

u_dt<-data.table(svd.out$u)
u_dt_melt<-melt(u_dt[,1:5])
ggplot(u_dt_melt)+geom_violin(aes(x=as.factor(variable),y=value),draw_quantiles=c(0.25,0.5,0.75))
#homework no 4 !!!!

#HHHHHHHHHHHHHHHHHHHHHHHHH

#HHHHHHHHHHHHHHHHHHHHHHHHH

k_means<-kmeans(edata,5)
clust<-k_means$cluster

library(irlba)
library(Rtsne)

# Set a seed for reproducible results
set.seed(1)
# complexity is a hyperparameter needed for this algorithm. 30 is a default
tsne_out <- Rtsne(edata,pca=TRUE,perplexity=30)
tsne_out_k = data.table(tsne_out$Y, clust)
ggplot(tsne_out_k) + geom_point(aes(x=V1, y=V2, col=as.factor(clust)))

#kmeans

#HHHHHHHHHHHHHHHHHHHHHHHHH


#role of Normalization

pc1<-prcomp(edata)
plot(pc1$rotation[,1],svd.out$v[,1])

edata.col<-scale(edata,scale=FALSE,center=TRUE)
svd.col<-svd(edata.col)
plot(pc1$rotation[,1],svd.col$v[,1],col=2)
abline(0,1)

library(irlba)
library(ggplot2)

tsvd.out<-irlba(edata,nv=4)
#nv is number of singular vectors
dim(tsvd.out$u)

plot(tsvd.out$v[,1],-1*svd.out$v[,1]); abline(0,1,col='red')

install.packages('Rtsne')
library(Rtsne)

set.seed(1)
tsne_out<-Rtsne(edata,pca=TRUE, perplexity=30)
tsne_out<-data.table(tsne_out$Y)
ggplot(tsne_out) + geom_point(aes(x=V1,y=V2))

#tSNE moze byc lepsze niz PCA
#w niektorych przypadkach

#homework 5 !!!
#use Kmeans, with k=5









