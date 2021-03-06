if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("Biobase")

library(Biobase)
library(GEOquery)

dat <- getGEO('GDS39', destdir = ".")

install.packages("data.table")
library(data.table)

geneexp <- Table(dat)
geneexp.tidy <- gather(geneexp, key="Samples", value="GeneExp", -c(1,2))

install.packages('reshape2')
library(reshape2)

dat.geneexp.complete_melted <- melt(dat.geneexp.complete)

homework_heatmap <- ggplot(data = dat.geneexp.complete_melted, mapping = aes(x = Var2, y = Var1, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = 'blue1', high = 'yellow1', midpoint = 0, limit = c(-3,3), space = "Lab")+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill="white"))+
  theme_void()

print(homework_heatmap)




