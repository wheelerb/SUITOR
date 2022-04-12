
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
    n.cores     \tab Number of cores to use for parallel computing \tab 1 \cr
    n.starts    \tab Number of starting points     \tab 30            \cr
    print       \tab 0 or 1 to print info         \tab 1             \cr 
    type        \tab Socket type in \code{\link[parallel]{makeCluster}} for 
        parallel computing \tab NULL       
    }

    \bold{Parallel computing} \cr
    If \code{n.cores > 1} and \code{type = NULL}, then \code{type} will be 
    set to "FORK" when running on unix, and set to "PSOCK" otherwise. \cr
    \bold{NOTE:} The R_LIBS_USER environment variable may need to be set to 
    the path where the SUITOR R package was installed, especially if it 
    was not installed in the default location. Examples of setting the 
    R_LIBS_USER environment variable on unix and windows are below. \cr
    (unix) export R_LIBS_USER=/path/to/SUITOR/package \cr
    (windows) set R_LIBS_USER=/path/to/SUITOR/package \cr
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
