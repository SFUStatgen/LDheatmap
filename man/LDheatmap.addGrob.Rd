\name{LDheatmap.addGrob}
\alias{LDheatmap.addGrob}
\title{
Add a graphical object to an LDheatmap plot
}
\description{
Add a graphical object to an LDheatmap plot such that the x-axis corresponds
to the physical map on the heatmap.
}
\usage{
LDheatmap.addGrob(LDheatmap, grob, height = 0.2)
}
\arguments{
  \item{LDheatmap}{
An object of class LDheatmap.
}
  \item{grob}{
A graphical object of class\code{\link[grid]{grob}}.
}
  \item{height}{
The height of the viewport in which the grob will be placed.
}
}
\value{
An object of class LDheatmap given as an argument, with the \code{grob}
\code{LDheatmapGrob} modified to inclue the new child grob.
}
\author{
Sigal Blay <sblay@sfu.ca>
}
\seealso{
\code{\link{LDheatmap}}
}
\examples{
% Add an empty rectangle frame
data(GIMAP5.CEU) 
ll<-LDheatmap(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)
llplusgrob<-LDheatmap.addGrob(ll, grid::rectGrob())
}
\keyword{ hplot }
