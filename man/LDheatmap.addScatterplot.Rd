\name{LDheatmap.addScatterplot}
\alias{LDheatmap.addScatterplot}
\title{
Add a scatter plot to an LDheatmap object
}
\description{
Add a scatter plot to an LDheatmap object. 
The x axis is the map of genetic distances of the SNPs.
}
\usage{
LDheatmap.addScatterplot(LDheatmap, P, height = 0.2, ylab = NULL, ylim=NULL,
type = "points")
}
\arguments{
  \item{LDheatmap}{
An object of class LDheatmap.
}
  \item{P}{
A vector with the values to be plotted as the y axis.
}
  \item{height}{
The height of the plot.
}
  \item{ylab}{
The y axis label.
}
  \item{ylim}{
The y axis limits.
}
  \item{type}{
Plot type. Possible values are \code{"points"} (the default), \code{"lines"} or \code{"both"}.
}
}
\details{
The function creates an \code{"association"} \code{grob} and adds it to the \code{LDheatmap} object.
Then it pushes a \code{viewport} and draws the \code{LDheatmapGrob} onto it.
}
\value{
An object of class LDheatmap given as an argument, with the \code{grob} 
\code{LDheatmapGrob} modified to inclue the \code{"association"} child grob.
}
\author{
Sigal Blay <sblay@sfu.ca> and more
}
\seealso{
\code{\link{LDheatmap}}
}
\examples{
# Load the package's data set
data("CEUData")
# Produce an LDheatmap object
MyLDheatmap <- LDheatmap(CEUSNP, genetic.distances = CEUDist, flip = TRUE)
# Generate an arbitrary vector of values to plot
Yvalues <- seq(length = length(MyLDheatmap$genetic.distances), from = 0.01, to = 0.5)
# Add scatter plot
assoc <- LDheatmap.addScatterplot(MyLDheatmap, Yvalues)
}
\keyword{ hplot }

