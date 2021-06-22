#setwd('/Users/jakub/Desktop/CBS_project/ESPCA/')

source('/Users/jakub/Desktop/CBS_project/ESPCA/fun_ESPCA.R')
source('/Users/jakub/Desktop/CBS_project/ESPCA/fun_SPCA.R')

BiocManager::install(c("fgsea"))
library(fgsea)
BiocManager::install(c("rstan"))
library(rstan)
BiocManager::install(c("ShortRead", 'GEOquery'))
library(ShortRead)
library(GEOquery)
# Loading ranks

ranks <- read.table('/Users/jakub/Desktop/CBS_project/ESPCA/brca_hd_tep_ranks.rnk', sep="\t", blank.lines.skip = TRUE, header = TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$rank, ranks$gene)
str(ranks)

# Loading pathways
pathways <- gmtPathways("/Users/jakub/Desktop/CBS_project/ESPCA/Human_GOBP_AllPathways_no_GO_iea_February_01_2017_symbol.gmt")
str(head(pathways))

# Running fgsea
gseaRes <- fgsea(pathways, ranks, minSize = 15, maxSize = 500)
edgess <- gseaRes$leadingEdge

# We need a graph with integers as nodes, not gene names
int_edgess = list()

for(i in edgess){
  for(j in i){
    int_edgess <- c(int_edgess, j)
  }
}

gene_names <- unique(int_edgess)
idx <- 1:length(gene_names)
names(idx) <- gene_names

the_edgess = edgess
counter = 1

for(i in edgess){
  new_names = c()
  for(j in i){
    new_name <- idx[j]
    new_names = c(new_names, new_name)
  }
  the_edgess[[counter]] = new_names
  counter = counter + 1
}

# Loading expr file
seqdata <- read.delim("/Users/jakub/Desktop/CBS_project/ESPCA/brca_hd_tep_expression.txt", stringsAsFactors = FALSE)
head(seqdata)

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
head(countdata)

#######
countdata.as.matrix <- as.matrix(countdata)
#######

# PCA
out11 = svd(countdata.as.matrix, nu = 2, nv = 2)

# SPCA
out22 =  SPCA(countdata.as.matrix, k = 2, kv = c(5,5))

# ESPCA
out33 = ESPCA(countdata.as.matrix, k = 2, the_edgess, k.group=2)

# Result 
PC.dat = cbind(cbind(out11$v, out22$V), out33$V)
colnames(PC.dat) = c("PCA.PC1","PCA.PC2","SPCA.PC1","SPCA.PC2","ESPCA.PC1","ESPCA.PC2")
row.names(PC.dat) = paste("Var",1:32, sep = "")
#++++++++++++++++++++++++++++++++++

PC1 <- data.table(out11$v)
PC2 <- data.table(out22$V)
PC3 <- data.table(out33$V)

pData = list()
pData[[1]] = names(countdata)
pData[[2]] = c('HD','HD','HD','HD','HD','HD','HD','HD','HD','HD','HD','HD','HD','HD','HD','HD','HD', 'BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA','BrCA')

PC1 <- data.table(out11$v, pData[[2]])
setnames(PC1, "V2", "donor")
setnames(PC1, "V2", "donor")
setnames(PC1, "donor", "V2")
PC1
PC2 <- data.table(out22$V, pData[[2]])
setnames(PC2, "V2", "donor")
setnames(PC2, "V2", "donor")
setnames(PC2, "donor", "V2")
PC3 <- data.table(out33$V, pData[[2]])
setnames(PC3, "V2", "donor")
setnames(PC3, "V2", "donor")
setnames(PC3, "donor", "V2")

ggplot(PC1)+geom_point(aes(x = V1, y = V2, col = as.factor(donor)))
ggplot(PC2)+geom_point(aes(x=V1,y=V2, col = as.factor(donor)))
ggplot(PC3)+geom_point(aes(x=V1,y=V2, col = as.factor(donor)))

#### Perc. of variance explained
# Normal PCA
U1 <- out11$u
v1 <- out11$v
D1 <- diag(out11$d)

plot(out11$d, type = 'o', )
abline(h=0, col = 'purple')

plot(out11$d^2/sum(out11$d^2)*100, type = 'o',ylab='Percent variability explained')
abline(h=0, col = 'purple')

d1 <- svd(countdata.as.matrix, nu = 2, nv = 2)$d
d1[1]^2/sum(d1^2)

var_explained1 <- d1
var_explained1 = var_explained1^2 / sum(var_explained1^2)
var_explained1 = as.data.frame(var_explained1)
var_explained1_1 <- var_explained1[1,1]
var_explained1_2 <- var_explained1[2,1]

ggplot(data = as.data.frame(var_explained1[1:2,]), aes(x = c(1:2), y = c(var_explained1_1, var_explained1_2))) + geom_bar(stat = 'identity', fill = 'steelblue') + theme_minimal() + xlab('Principal Component') + ylab('Perc. variance explained') + ylim(0,1) + ggtitle('Normal PCA')

# SPCA
U2 <- out22$U
v2 <- out22$V
D2 <- out22$D

d2 <- diag(out22$D)
d2[1]^2/sum(d2^2)

var_explained2 <- diag(D2)
var_explained2 = var_explained2^2 / sum(var_explained2^2)
var_explained2 = as.data.frame(var_explained2)
var_explained2_1 <- var_explained2[1,1]
var_explained2_2 <- var_explained2[2,1]

ggplot(data = as.data.frame(var_explained2[1:2,]), aes(x = c(1:2), y = c(var_explained2_1, var_explained2_2))) + geom_bar(stat = 'identity', fill = 'steelblue') + theme_minimal() + xlab('Principal Component') + ylab('Perc. variance explained') + ylim(0,1) + ggtitle('SPCA')

# ESPCA
U3 <- out33$U
v3 <- out33$V
D3 <- out33$D

d3 <- out33$D
d3[1]^2/sum(d3^2)

var_explained3 <- diag(D3)
var_explained3 = var_explained3^2 / sum(var_explained3^2)
var_explained3 = as.data.frame(var_explained3)
var_explained3_1 <- var_explained3[1,1]
var_explained3_2 <- var_explained3[2,1]

ggplot(data = as.data.frame(var_explained3[1:2,]), aes(x = c(1:2), y = c(var_explained3_1, var_explained3_2))) + geom_bar(stat = 'identity', fill = 'steelblue') + theme_minimal() + xlab('Principal Component') + ylab('Perc. variance explained') + ylim(0,1) + ggtitle('ESPCA')

# All plots in one (all values for 2 PC from all 3 techniques)
library(reshape2)
PC.dat.melt <- melt(PC.dat)

ggplot(data = PC.dat.melt, aes(x = Var2, y = value, fill = Var2)) + geom_bar(stat = 'identity', position = position_dodge()) + theme_minimal() +ylab('Value') + xlab('Principal Components')

# All variance explained plots in one plot
all_var_explained <- data.frame(x = c('PCA.PC1', 'PCA.PC2', 'SPCA.PC1', 'SPCA.PC2', 'ESPCA.PC1', 'ESPCA.PC2'), y = c(var_explained1_1, var_explained1_2, var_explained2_1, var_explained2_2, var_explained3_1, var_explained3_2))

ggplot(data = all_var_explained, aes(x = x, y = y, fill = x)) + geom_bar(stat = 'identity', position = position_dodge()) + theme_minimal() + xlab('Principal Components') + ylab('Variance explained') + ggtitle('PCA vs SPCA vs ESPCA comparison given the variance explained') + geom_hline(yintercept = 0.97261430) + geom_hline(yintercept = 0.02738570)









