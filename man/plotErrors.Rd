
\name{plotErrors}
\alias{plotErrors}
\title{Plot train and test errors}
\description{
    Plot train and test errors                       
}
\usage{ plotErrors(x)
}
\arguments{
    \item{x}{Data frame of \code{summary} results in the return list from 
        \code{\link{suitor}} or from \code{\link{getSummary}}, or
        a data frame with columns \code{Rank}, \code{Type}, and 
        \code{MSErr}.}
}
\details{
    The optimal rank is the minimum at which the test error is attained, 
    and appears as a red dot on the graph.
}
\value{ NULL }
%\references{}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
    data(plotData, package="SUITOR")
    plotErrors(plotData)
}
\keyword{ NMF, cross-validation, mutational signatures } 
%\seealso{ \code{\link{getSummary}} }
