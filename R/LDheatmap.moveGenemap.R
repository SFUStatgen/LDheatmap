# Genemap moving function
  # Test movement with symbols on

LDheatmap.moveGenemap <- function(LDheatmap, distance = "close"){
  if(is.null(LDheatmap$flipVP)){
    if(distance == "close") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.52, "snpc"), y = unit(0.48, "snpc")))
    if(distance == "medium") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.56, "snpc"), y = unit(0.44, "snpc")))
    if(distance == "far") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.60, "snpc"), y = unit(0.4, "snpc")))
  }
  else{
    if(distance == "close") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.5, "snpc"), y = unit(0.52, "snpc")))
    if(distance == "medium") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.5, "snpc"), y = unit(0.56, "snpc")))
    if(distance == "far") temp <- editGrob(LDheatmap$LDheatmapGrob, gPath("geneMap"), vp = viewport(x = unit(0.5, "snpc"), y = unit(0.6, "snpc")))
  }
  LDheatmap$LDheatmapGrob <- temp
  grid.newpage()
  grid.draw(LDheatmap$LDheatmapGrob)
  return(LDheatmap)
}
