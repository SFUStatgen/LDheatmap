library(LDheatmap)

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("snpStats")
biocLite("rtracklayer")

library(snpStats)
library(rtracklayer)
library(grid)


########################### Sample Code from Vignette ##########################
data(GIMAP5.CEU)
load(system.file("extdata/addTracks.RData",package="LDheatmap"))

ll <-LDheatmap(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)
# addGenes sometimes throws "Error: Bad Request" when called with the same parameters as addRecombRate, even though it grabs the requested values
# Just run grid.draw(llGenes$LDheatmapGrob) after to get the desired chart. Only necessary on my laptop for some reason
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
####################################################################

########################## RESOLVED ##############################
# Bug report: "Symbols" are gone when LDheatmap is generated with flip = TRUE
# Fix: LDheatmapMap.add() had a code section only if(flip) that set symbols <- NULL. 
    # Resolved by using same assignment as was used in the non-flipped case
require('LDheatmap')

 # taken from examples LDheatmap\demo\LDheatmap.R, works correctly
data(CEUData)
MyHeatmap <- LDheatmap(CEUSNP, genetic.distances = CEUDist,
                          color = grey.colors(20))
LDheatmap(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"))
childNames(grid.get("geneMap"))


 #BUG: With flip=TRUE it does not work anymore
  MyHeatmap <- LDTest(CEUSNP, genetic.distances = CEUDist,
                            color = grey.colors(20), flip=TRUE)
  MyHeatmapTest <- LDTest(CEUSNP, genetic.distances = CEUDist,
                          color = grey.colors(20), flip=FALSE)
  
# Previous code  
 LDheatmap(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"))
 childNames(grid.get("geneMap")) # No symbols
 LDheatmap(MyHeatmapTest, SNP.name = c("rs2283092", "rs6979287"))
 childNames(grid.get("geneMap")) # Symbols available only in the non-flipped case
# New code
 LDTest(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"), text = FALSE)
 
 childNames(grid.get("geneMap")) # Symbols available
 LDTest(MyHeatmapTest, SNP.name = c("rs2283092", "rs6979287"))
 childNames(grid.get("geneMap")) # Symbols available
##################################################################

####################### UNRESOLVED ###########################
# Bug report: Cannot change the heatmap position away from snp names and segment. 
 # Other questions: Is there a way to add bp positions? (i.e. beginning and end positions at the segment bar)
 # When using grid.edit(gPath("LDheatmap", "heatMap", "heatmap"), gp = gpar(col = "gray90", lwd = 1)) with a colour other than white,
 # the whole square shows up. How can this be done only on the "heatmap" side of the plot?
testMap <- LDTest(CEUSNP, genetic.distances = CEUDist,
                               color = grey.colors(20), flip = TRUE)
LDTest(testMap, SNP.name = c("rs2283092", "rs6979287"))
grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=1.5))
new.grob <- editGrob(testMap$LDheatmapGrob, gPath("geneMap", "segments"),
                     gp=gpar(col="orange"))
moveGenemap <- editGrob(testMap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.55, "snpc"), y = unit(0.45, "snpc")))
# We can control how far from the heatmap the geneMap is by fixing the vp above with x = unit(0.5, "snpc") and the y value 
# as something larger than 0.49
grid.newpage()
grid.draw(moveGenemap)

LDTest.moveGenemap <- function(LDheatmap, distance = "close"){
  if(is.null(LDheatmap$flipVP)){
    print("Not Flipped")
    if(distance == "close") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.52, "snpc"), y = unit(0.48, "snpc")))
    if(distance == "medium") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.56, "snpc"), y = unit(0.44, "snpc")))
    if(distance == "far") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.60, "snpc"), y = unit(0.4, "snpc")))
  }
  else{
    print("Flipped")
    if(distance == "close") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.5, "snpc"), y = unit(0.52, "snpc")))
    if(distance == "medium") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.5, "snpc"), y = unit(0.56, "snpc")))
    if(distance == "far") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.5, "snpc"), y = unit(0.6, "snpc")))
  }
  LDheatmap$LDheatmapGrob <- temp
  return(LDheatmap)
}

movedGenemap <- LDTest.moveGenemap(testMap, distance = "medium")
grid.newpage()
grid.draw(movedGenemap$LDheatmapGrob)
##############################################################

############################### RESOLVED ####################################
# Bug report: highlight function needs to be mirrored when used with flip = TRUE #
# Sample of highlight on normal heatmap
library(mvtnorm)
data(CEUData)
tt <- LDheatmap(CEUSNP, genetic.distances=CEUDist)
LDheatmap.highlight(tt, 3, 8, col="blue", fill="green", lwd=3)

# Sample of highlight on flipped heatmap
ttFlip <- LDheatmap(CEUSNP, genetic.distances = CEUDist, flip = TRUE)
LDheatmap.highlight(ttFlip, 6, 9, col = "blue", fill = "green", lwd = 3)

# Solution: Add flipOutline parameter such that if a user encounters a flip problem they can change the flipOutline param value to reverse
LDheatmap.highlight(ttFlip, 6, 9, col = "blue", fill = "green", lwd = 3, flipOutline = TRUE)
#############################################################################



################### RESOLVED ####################################
# Bug report: Text issue where the text does not follow the heatmap during flip
# Fix: Added flip parameter to makeImageText() such that the data could be added to the graph in the appropriate order. Adjusted the call
    # of makeImageText() in LDheatmap() to accommodate for the added parameter.
llText <- LDTest(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE, text = TRUE)

###################################################################


################### UNRESOLVED ####################################
# Bug report: Using 1 or 2 scatterplots above the heatmap work but the third one causes problems, doesnt display appropriately
library(LDheatmap)
library(viridis)
data(GIMAP5.CEU)
load(system.file("extdata/addTracks.RData",package="LDheatmap"))

ll2 <- LDheatmap(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE,color = viridis(20))
llGenes2 <- LDheatmap.addGenes(ll, chr="chr7", genome="hg18")

grid.newpage()
grid.draw(llGenes2$LDheatmapGrob)
llGenesRecomb <- LDheatmap.addRecombRate(llGenes, chr="chr7", genome="hg18")
grid.newpage()
grid.draw(llGenesRecomb$LDheatmapGrob)

set.seed(1)
atests<-runif(nrow(GIMAP5.CEU$snp.support))
names(atests)<-rownames(GIMAP5.CEU$snp.support)
atests["rs6598"]<-1e-5

llGenesRecombScatter<-LDheatmap.addScatterplotztw(llGenes,-log10(atests), ylab="GWAS\n-log10(p-values)")
llGenesRecombScatter2<-LDheatmap.addScatterplotztw2(llGenesRecombScatter,-log10(atests), ylab="eQTL \n-log10(p-values)")
#pdf('tmp.pdf',width = 10,height = 8)
llGenesRecombScatter3<-LDheatmap.addScatterplotztw3(llGenesRecombScatter2,-log10(atests), ylab="CLPP")

grid.newpage()
pushViewport(llGenesRecombScatter3$flipVP)
grid.rect(gp = gpar(col = "red"))

upViewport()
pushViewport(llGenesRecombScatter3$heatmapVP)
grid.rect(gp = gpar(col = "blue"))

upViewport()
pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association$vp)
grid.rect(gp = gpar(col = "black"))

# Not presented properly because they have their vp nested
upViewport()
pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association2$vp)
grid.rect(gp = gpar(col = "green"))

# Not presented properly because they have their vp nested
upViewport()
pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association3$vp)
grid.rect(gp = gpar(col = "pink"))


grid.newpage()
pushViewport(llGenesRecombScatter3$flipVP)
grid.rect(gp = gpar(col = "red"))

pushViewport(llGenesRecombScatter3$LDheatmapGrob$vp)
grid.rect(gp = gpar(col = "blue"))

#pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association$vp)
grid.rect(gp = gpar(col = "black"))
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$association)

#pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association2$vp)
pushViewport(llGenesRecombScatter3$flipVP)
grid.rect(gp = gpar(col = "red"))

pushViewport(llGenesRecombScatter3$LDheatmapGrob$vp)
grid.rect(gp = gpar(col = "blue"))
grid.rect(gp = gpar(col = "green"))
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$association2)

#pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association3$vp)
grid.rect(gp = gpar(col = "pink"))
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$association3)


vp <- constructVP(llGenes$LDheatmapGrob, 0.03, TRUE)
vpStacked <- vpStack(llGenes$flipVP,vp)

grid.newpage()
pushViewport(llGenes$flipVP)
pushViewport(vp)
grid.rect(gp = gpar(col = "red"))
upViewport()
pushViewport(vpStacked)
grid.rect(gp = gpar(col = "blue"))

grid.ls(llGenes$LDheatmapGrob, grobs = FALSE, viewports = TRUE)
dev.off()

grid.newpage()
pushViewport(llGenesRecombScatter3$heatmapVP)
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$heatMap)
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$geneMap)
pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$transcripts$vp)
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$transcripts$children)
upViewport()
upViewport()
pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association$vp)
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$association$children)
upViewport()
upViewport()
pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association2$vp)
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$association2$children)
upViewport()
upViewport()
pushViewport(llGenesRecombScatter3$LDheatmapGrob$children$association3$vp)
grid.draw(llGenesRecombScatter3$LDheatmapGrob$children$association3$children)

grid.newpage()
grid.draw(llGenesRecombScatter2$LDheatmapGrob)

# 
# posn<-GIMAP5.CEU$snp.support$Position
# manhattan2<-ggplotGrob( {
#   qplot(posn,-log10(atests),xlab="", xlim=range(posn),asp=1/10)
#   last_plot() + theme(axis.text.x=element_blank(), axis.title.y = element_text(size = rel(0.75)))
# }
# )
# pdf('tmp.pdf',width = 10,height = 8)
# llQplot<-LDheatmap.addGrob(ll,manhattan2,)
# dev.off()
# 
# manhattan2<-ggplotGrob(
#   {
#     qplot(posn,-log10(atests),xlab="", xlim=range(posn),asp=1/10)
#     last_plot() + theme(axis.text.x=element_blank(),axis.title.y = element_text(size = rel(0.75)))
#   }
# )
# 
# pdf('tmp.pdf',width = 10,height = 8)
# llQplot2<-LDheatmap.addGrob(llGenesRecombScatter,rectGrob(gp=gpar(col="white")),height=.2)
# pushViewport(viewport(x=.485,y=.775,width=.775,height=.9))
# grid.draw(manhattan2)
# popViewport(1)
# dev.off()
# 
# 



##### define color #####

LDheatmap.addScatterplotztw <- function (LDheatmap, P, height = 0.2, ylab = NULL, ylim = NULL,  type = "points",color="red",pch=16) 
{
  if (dim(LDheatmap$LDmatrix)[1] != length(P)) {
    print("Length of vector not equal number of SNPs in LDheatmap")
    return()
  }
  flip <- !is.null(LDheatmap$flipVP)
  vp <- constructVP(LDheatmap$LDheatmapGrob, 0.03, flip)
  vp$height <- unit(height, "npc")
  vp$name <- "associationVP"
  if (is.null(ylim)) 
    ylim <- c(floor(min(P)), ceiling(max(P)))
  vp$yscale <- ylim
  vp$xscale <- c(min(LDheatmap$genetic.distances), max(LDheatmap$genetic.distances))
  xaxis <- linesGrob(x = vp$xscale, y = 0, default.units = "native", 
                     name = "xaxis")
  yaxis <- linesGrob(x = min(LDheatmap$genetic.distances), 
                     y = vp$yscale, default.units = "native", name = "yaxis")
  yaxisT <- yaxisGrob(name = "yaxis_ticks", gp = gpar(fontsize = 7))
  ylab <- textGrob(ylab, rot = 90, gp = gpar(fontsize = 9), 
                   name = "yaxis_title", x = unit(min(LDheatmap$genetic.distances), 
                                                  "native") - unit(10, "millimeters"))
  vpstack <- vp
  if (flip) 
    vpstack <- vpStack(LDheatmap$flipVP, vp)
  association <- gTree(children = gList(xaxis, yaxis, yaxisT, 
                                        ylab), name = "association", vp = vpstack)
  if (type == "points" || type == "both") {
    graph_points <- pointsGrob(LDheatmap$genetic.distances, P, size = unit(2, "millimeters"), name = "points",pch=pch, gp=gpar(col=color))
    association <- addGrob(association, graph_points)
  }
  if (type == "lines" || type == "both") {
    graph_lines <- linesGrob(LDheatmap$genetic.distances, 
                             P, default.units = "native", name = "lines")
    association <- addGrob(association, graph_lines)
  }
  LDheatmap$LDheatmapGrob <- addGrob(LDheatmap$LDheatmapGrob, 
                                     association)
  LDheatmap$LDheatmapGrob <- moveTitles(LDheatmap$LDheatmapGrob, 
                                        vp)
  return(LDheatmap)
}

environment(LDheatmap.addScatterplotztw) <- asNamespace('LDheatmap')




##### define color2 ####
LDheatmap.addScatterplotztw2 <- function (LDheatmap, P, height = 0.2, ylab = NULL, ylim = NULL,  type = "points",color="purple",pch=16) 
{
  if (dim(LDheatmap$LDmatrix)[1] != length(P)) {
    print("Length of vector not equal number of SNPs in LDheatmap")
    return()
  }
  flip <- !is.null(LDheatmap$flipVP)
  vp <- constructVP(LDheatmap$LDheatmapGrob, 0.03, flip)
  vp$height <- unit(height, "npc")
  vp$name <- "associationVP"
  if (is.null(ylim)) 
    ylim <- c(floor(min(P)), ceiling(max(P)))
  vp$yscale <- ylim
  vp$xscale <- c(min(LDheatmap$genetic.distances), max(LDheatmap$genetic.distances))
  xaxis <- linesGrob(x = vp$xscale, y = 0, default.units = "native", 
                     name = "xaxis")
  yaxis <- linesGrob(x = min(LDheatmap$genetic.distances), 
                     y = vp$yscale, default.units = "native", name = "yaxis")
  yaxisT <- yaxisGrob(name = "yaxis_ticks", gp = gpar(fontsize = 7))
  ylab <- textGrob(ylab, rot = 90, gp = gpar(fontsize = 9), 
                   name = "yaxis_title", x = unit(min(LDheatmap$genetic.distances), 
                                                  "native") - unit(10, "millimeters"))
  vpstack <- vp
  if (flip) 
    vpstack <- vpStack(LDheatmap$flipVP, vp)
  association2 <- gTree(children = gList(xaxis, yaxis, yaxisT, 
                                         ylab), name = "association2", vp = vpstack)
  if (type == "points" || type == "both") {
    graph_points <- pointsGrob(LDheatmap$genetic.distances, P, size = unit(2, "millimeters"), name = "points",pch=pch, gp=gpar(col=color))
    association2 <- addGrob(association2, graph_points)
  }
  if (type == "lines" || type == "both") {
    graph_lines <- linesGrob(LDheatmap$genetic.distances, 
                             P, default.units = "native", name = "lines")
    association2 <- addGrob(association2, graph_lines)
  }
  LDheatmap$LDheatmapGrob <- addGrob(LDheatmap$LDheatmapGrob, 
                                     association2)
  LDheatmap$LDheatmapGrob <- moveTitles(LDheatmap$LDheatmapGrob, 
                                        vp)
  return(LDheatmap)
}

environment(LDheatmap.addScatterplotztw2) <- asNamespace('LDheatmap')





##### define color3 ####
LDheatmap.addScatterplotztw3 <- function (LDheatmap, P, height = 0.2, ylab = NULL, ylim = NULL,  type = "points",color="green4",pch=16) 
{
  if (dim(LDheatmap$LDmatrix)[1] != length(P)) {
    print("Length of vector not equal number of SNPs in LDheatmap")
    return()
  }
  flip <- !is.null(LDheatmap$flipVP)
  vp <- constructVP(LDheatmap$LDheatmapGrob, 0.03, flip)
  vp$height <- unit(height, "npc")
  vp$name <- "associationVP"
  if (is.null(ylim)) 
    ylim <- c(floor(min(P)), ceiling(max(P)))
  vp$yscale <- ylim
  vp$xscale <- c(min(LDheatmap$genetic.distances), max(LDheatmap$genetic.distances))
  xaxis <- linesGrob(x = vp$xscale, y = 0, default.units = "native", 
                     name = "xaxis")
  yaxis <- linesGrob(x = min(LDheatmap$genetic.distances), 
                     y = vp$yscale, default.units = "native", name = "yaxis")
  yaxisT <- yaxisGrob(name = "yaxis_ticks", gp = gpar(fontsize = 7))
  ylab <- textGrob(ylab, rot = 90, gp = gpar(fontsize = 9), 
                   name = "yaxis_title", x = unit(min(LDheatmap$genetic.distances), 
                                                  "native") - unit(10, "millimeters"))
  vpstack <- vp
  if (flip) 
    vpstack <- vpStack(LDheatmap$flipVP, vp)
  association3 <- gTree(children = gList(xaxis, yaxis, yaxisT, 
                                         ylab), name = "association3", vp = vpstack)
  if (type == "points" || type == "both") {
    graph_points <- pointsGrob(LDheatmap$genetic.distances, P, size = unit(2, "millimeters"), name = "points",pch=pch, gp=gpar(col=color))
    association3 <- addGrob(association3, graph_points)
  }
  if (type == "lines" || type == "both") {
    graph_lines <- linesGrob(LDheatmap$genetic.distances, 
                             P, default.units = "native", name = "lines")
    association3 <- addGrob(association3, graph_lines)
  }
  LDheatmap$LDheatmapGrob <- addGrob(LDheatmap$LDheatmapGrob, 
                                     association3)
  LDheatmap$LDheatmapGrob <- moveTitles(LDheatmap$LDheatmapGrob, 
                                        vp)
  return(LDheatmap)
}

environment(LDheatmap.addScatterplotztw3) <- asNamespace('LDheatmap')
