library(LDheatmap)

source("https://bioconductor.org/biocLite.R")
biocLite("snpStats")
biocLite("rtracklayer")

library(snpStats)
library(rtracklayer)
library(grid)

# Bug report testing
require('LDheatmap')

 # taken from examples LDheatmap\demo\LDheatmap.R, works correctly
data(CEUData)
MyHeatmap <- LDheatmap(CEUSNP, genetic.distances = CEUDist,
                          color = grey.colors(20))
LDheatmap(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"))
childNames(grid.get("geneMap"))


 #BUG: With flip=TRUE it does not work anymore
  MyHeatmap <- LDheatmap(CEUSNP, genetic.distances = CEUDist,
                            color = grey.colors(20), flip=FALSE)
 LDheatmap(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"))
 childNames(grid.get("geneMap"))



data(GIMAP5.CEU)
load(system.file("extdata/addTracks.RData",package="LDheatmap"))

ll <-LDheatmap(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)
llGenes <- LDheatmap.addGenes(ll, chr="chr7", genome="hg18")

grid.newpage()
grid.draw(llGenes$LDheatmapGrob)

llGenesRecomb <- LDheatmap.addRecombRate(llGenes, chr="chr7", genome="hg18")

grid.newpage()
grid.draw(llGenesRecomb$LDheatmapGrob)

set.seed(1)
atests <- runif(nrow(GIMAP5.CEU$snp.support))
names(atests) <- rownames(GIMAP5.CEU$snp.support)
atests["rs6598"]<-1e-5

llGenesRecombScatter <- LDheatmap.addScatterplot(llGenesRecomb,-log10(atests),
                                                 ylab="-log10(p-values)")

# Adds Manhatten plot above our graph using qplot function in ggplot
require("ggplot2")
posn <- GIMAP5.CEU$snp.support$Position
manhattan2<-ggplotGrob(
     {
      qplot(posn,-log10(atests),xlab="", xlim=range(posn),asp=1/10)
      last_plot() + theme(axis.text.x=element_blank(),
                           axis.title.y = element_text(size = rel(0.75)))
     }
   )
llQplot<-LDheatmap.addGrob(ll,manhattan2,height=.7)

llImage<-LDheatmap.addGrob(ll,rasterGrob(GIMAP5ideo))
