# LDheatmap testing functionality


# source("https://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("snpStats")


library(LDheatmap)
library(grid)
library(snpStats)

data(GIMAP5.CEU)
ll<-testFN(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)

dim(GIMAP5.CEU$snp.data)[2] # Number of rows
colcutMatrixByCol <- matrix(globalColcut, ncol = 23, byrow = FALSE)
colcutMatrixByCol
colcutMatrixByRow <- matrix(globalColcut, nrow = 23, byrow = TRUE)
colcutMatrixByRow
# It appears that filling byrow versus bycol actually swaps from upper triangular to lower triangular matrix
# This makes sense as the last column is the same as the last row when using the opposite methods b/c square matrix
grid.newpage()
grid.draw(globalRect) # When flip = TRUE, this draws under the diagonal from bottom left to top right in the grid page
                      # When flip = FALSE this draws above the diagonal from bottom left to top right in the grid page
grid.newpage()
grid.draw(globalRotatedRect)

# Testing alternative ways
dimension <- 23
grid.lines()
grid.newpage()

# Diagonals to the right
xcoords <- (1:dimension)/dimension
xend <- rep(1, dimension)
ycoords <- ((dimension-1):0)/dimension
ystart <- rep(0, dimension)
testX <- c(rbind(xcoords, xend))
testY <- c(rbind(ystart, ycoords))

# Diagonals to the left
xcoords2 <- (dimension:1)/dimension
xend2 <- rep(0, dimension)
ycoordsLeft <- (dimension:1)/dimension
testX2 <- c(rbind(xcoords2, xend2))
testY2 <- c(rbind(ystart, ycoordsLeft))

grid.polyline(x = c(testX, testX2),
              y = c(testY, testY2),
              id.lengths = rep(2, dimension*2),
              gp = gpar(col =1)
)

# Image is upside down, try to replace all with 1- val
flipTestX <- 1- testX
flipTestX2<- 1- testX2
flipTestY <- 1- testY
flipTestY2<- 1- testY2
grid.newpage()
grid.polyline(x = c(flipTestX, flipTestX2),
              y = c(flipTestY, flipTestY2),
              id.lengths = rep(2, dimension*2),
              gp = gpar(col = 1)
              )

# Try this again but use 0.5, seeing if we can shift the now flipped image downwards
flipTestHalfY <- 0.5- testY
flipTestHalfY2<- 0.5- testY2
grid.newpage()
grid.polyline(x = c(flipTestX, flipTestX2),
              y = c(flipTestHalfY, flipTestHalfY2),
              id.lengths = rep(2, dimension*2),
              gp = gpar(col = 1)
)

# Image is pretty much where we want it, now invesigate to see if we can colour segments


# Other idea: Subset one vp (angled) within the other (not-angled) and see if tracklaying is easier
grid.newpage()
vp <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"), x = 0.5, y = 0.4,
                name="base")
pushViewport(vp)
grid.rect(gp=gpar(lty="dashed", col = "black"))
rotatedVP <- viewport(width = unit(0.8, "snpc"), height= unit(0.8, "snpc"), x = 0.5, y = 0.5, angle=-45, name="rotate", clip = FALSE)
pushViewport(rotatedVP)
grid.rect(gp=gpar(lty="dashed", col = "red")) #, just = c("left", "top"))
upViewport()
subsetVP <- viewport(width = unit(1.13, "snpc"), height = unit(0.4, "snpc"), x = 0.5, y = 0.8)
pushViewport(subsetVP)
grid.rect(gp=gpar(lty="dashed", col = "purple"))
trackVP <- viewport(width = unit(0.8, "snpc"), height = unit(0.4, "snpc"), x= 0.5, y =0.75)

# Try to add these viewports to the graphical output of ldheatmap
grid.newpage()
pushViewport(vp)
grid.draw(ll$LDheatmapGrob)
grid.rect(gp=gpar(lty="dashed", col = "black"))
upViewport()
pushViewport(trackVP)
grid.rect(gp=gpar(lty = "dashed", col = "purple"))

# Same thing but with track layer
# Currently not working, review
library(rtracklayer)
llplusgenes <- LDheatmap.addGenes(ll, chr="chr7", genome="hg18")
grid.newpage()
pushViewport(vp)
grid.draw(llplusgenes$LDheatmapGrob)
grid.rect(gp=gpar(lty="dashed", col = "black"))
upViewport()
pushViewport(trackVP)
grid.rect(gp=gpar(lty = "dashed", col = "purple"))



# What about just drawing heatmaps with ggplot instead?
library(reshape2)
library(ggplot2)
library(dplyr) # For restructuring

# Same data as colcut but as numeric rather than a cut factor
valueMatrixByCol <- matrix(globalImgLDmatrix, ncol = 23, byrow = FALSE) %>% as.data.frame() %>% mutate_all(as.numeric) 
valueMatrixByCol
valueMatrixByRow <- matrix(globalImgLDmatrix, nrow = 23, byrow = TRUE) %>% as.data.frame() %>% mutate_all(as.numeric)
valueMatrixByRow

# Generate DF structures, use character data so that no melt problems if feeding in colcut. If using numeric, no problem
colDF <- as.data.frame(cbind(names(GIMAP5.CEU$snp.data@.Data[1,]), valueMatrixByCol))
colDFChar <- as.data.frame(cbind(names(GIMAP5.CEU$snp.data@.Data[1,]), colcutMatrixByCol)) %>% mutate_all(as.character)
rowDF <- as.data.frame(cbind(rev(names(GIMAP5.CEU$snp.data@.Data[1,])), valueMatrixByRow))
#rowDF <- rowDF %>% mutate_all(as.character)

# Need column names
names(colDF) <- c("genes", names(GIMAP5.CEU$snp.data@.Data[1,]))
names(colDFChar) <- c("genes", names(GIMAP5.CEU$snp.data@.Data[1,]))
names(rowDF) <- c("genes", rev(names(GIMAP5.CEU$snp.data@.Data[1,])))

# Melt into long format, can convert back to factor if desired
meltedColDF <- melt(colDF, id.vars = "genes")
meltedColDFChar <- melt(colDFChar, id.vars = "genes")
meltedRowDF <- melt(rowDF, id.vars = "genes")

# Plot the heatmap
# Notes: Can change scale_fill_continuous to scale_fill_discrete and use the colour cut values
  # If problems occur changing the colour in this format, see: https://stackoverflow.com/questions/6906661/ggplot2-make-missing-value-in-geom-tile-not-blank
  # Currently has a gray border around plot, need to remove
library(forcats) # Used to order x and y axis accordingly

# This is the eventual flipped plot though it has not been rotated yet
flipButNotRotated <- ggplot(meltedColDF, aes(variable, fct_inorder(genes))) +
  geom_tile(aes(fill = value), color = "white", na.rm = TRUE) +
  scale_fill_continuous(na.value = "white", low = "yellow", high = "red") +
  ylab("List of genes ") +
  xlab("List of variables (genes)") +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.position = c(0.2, 0.8),
        plot.title = element_text(size=16),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        aspect.ratio = 1) +
  labs(fill = "Expression level") 

#Extract Legend 
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(flipButNotRotated)
flipButNotRotated <- flipButNotRotated + theme(legend.position ="none")

# Figure out if we can use push and popviewport to travel appropriately through ggplot objects instead of just using print
library(grid)
grid.newpage()
#testVP <- viewport(angle = -45, width = unit(0.7, "snpc"), height = unit(0.7, "snpc"), clip = FALSE, name = "test")
#grid.rect()
library(ggplotify)
# there is an alternative as.ggplot() that may be useful for tidying this up and testing in the other graphical format
holder <- as.grob(flipButNotRotated)
grid.newpage()
# The plot below has x axis top left to bottom right diagonal, y axis bottom left to top right
globalGenemap$children$diagonal$vp <- viewport(x = 0.45, y = 0.45, width = unit(0.8, "snpc"), height = unit(0.8, "snpc"), angle = 0)
globalGenemap$children$segments$vp <- viewport(x = 0.45, y = 0.45, width = unit(0.8, "snpc"), height = unit(0.8, "snpc"), angle = 0)
globalGenemap$children$title$rot <- 45
treeTest <- gTree(children = gList(holder, globalGenemap),
                  vp = testVP)
grid.draw(treeTest)

grid.draw(flipButNotRotated, vp = "test")
print(flipButNotRotated, vp = viewport(angle = -45, width = unit(0.7, "snpc"), height = unit(0.7, "snpc"), clip = FALSE))
grid.draw(legend)

# Try to add the genemap now
# genemapVP <- viewport(width = unit(0.8, "snpc"), height = unit(0.8, "snpc"), name = "genemapVP")
# globalGenemap$children$diagonal$vp$width <- unit(0.8, "snpc")
# globalGenemap$children$diagonal$vp$height <- unit(0.8, "snpc")
print(globalGenemap$children$segments, vp = viewport(width = unit(0.7, "snpc"), height = unit(0.7, "snpc")))

# This is the flip = false plot
noFlipNoRotate <- ggplot(meltedRowDF, aes(variable, fct_inorder(genes))) +
  geom_tile(aes(fill = value), color = "white", na.rm = TRUE) +
  scale_fill_continuous(na.value = "white", low = "yellow", high = "red") +
  ylab("List of genes ") +
  xlab("List of variables (genes)") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.position = c(0.8, 0.2),
        plot.title = element_text(size=16),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        aspect.ratio = 1) +
  labs(fill = "Expression level") 

noFlipNoRotate


clip.demo <- function(i, j, clip1, clip2) {
  pushViewport(viewport(layout.pos.col=i,
                        layout.pos.row=j))
  pushViewport(viewport(width=0.6, height=0.6, clip=clip1))
  grid.rect(gp=gpar(fill="white"))
  grid.circle(r=0.55, gp=gpar(col="red", fill="pink"))
  popViewport()
  pushViewport(viewport(width=0.6, height=0.6, clip=clip2))
  grid.polygon(x=c(0.5, 1.1, 0.6, 1.1, 0.5, -0.1, 0.4, -0.1),
               y=c(0.6, 1.1, 0.5, -0.1, 0.4, -0.1, 0.5, 1.1),
               gp=gpar(col="blue", fill="light blue"))
  popViewport(2)
}

grid.newpage()
grid.rect(gp=gpar(fill="grey"))
pushViewport(viewport(layout=grid.layout(2, 2)))
clip.demo(1, 1, FALSE, FALSE)
clip.demo(1, 2, TRUE, FALSE)
clip.demo(2, 1, FALSE, TRUE)
clip.demo(2, 2, TRUE, TRUE)
popViewport()

pushViewport(plotViewport(c(5, 4, 2, 2)))
pushViewport(dataViewport(pressure$temperature,
                            pressure$pressure,
                            name="plotRegion"))
grid.points(pressure$temperature, pressure$pressure,
            name="dataSymbols")
grid.rect()
grid.xaxis()
grid.yaxis()
