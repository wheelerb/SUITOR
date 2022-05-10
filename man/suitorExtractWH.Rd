
\name{suitorExtractWH}
\alias{suitorExtractWH}
\title{suitorExtractWH}
\description{
    Extract the matrix of activities (exposures) and matrix of signatures
}
\usage{ suitorExtractWH(data, rank, op=NULL)
}
\arguments{
    \item{data}{Data frame or matrix containing mutational signatures.
        This object must contain non-negative values}
    \item{rank}{Integer > 0}
    \item{op}{List of options (see details). The default is NULL.}
}
\details{
    \bold{Options list op:}
    \tabular{lll}{
    \bold{Name} \tab \bold{Description} \tab \bold{Default Value} \cr
    min.value   \tab Minimum value of matrix before factorizing \tab 1e-4 \cr
    BPPARAM    \tab See \code{\link[BiocParallel]{BiocParallelParam}} 
    \tab NULL \cr
    n.starts    \tab Number of starting points     \tab 30            \cr
    print       \tab 0 or 1 to print info         \tab 1             \cr       
    }

    \bold{Parallel computing} \cr
    The \code{BiocParallel} package is used for parallel computing.
    If \code{BPPARAM = NULL}, then \code{BPPARAM} 
    will be set to \code{\link[BiocParallel]{SerialParam}}. \cr

}
\value{ A list containing the objects:
    \itemize{
    \item \code{H}: Matrix of activities (exposures)
    \item \code{W}: Matrix of signatures
    }
}
%\references{}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
    data(SimData, package="SUITOR")

    suitorExtractWH(SimData, 2)
}
\keyword{ NMF, cross-validation, mutational signatures } 
\seealso{ \code{\link{suitor}} }
