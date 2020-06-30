
\name{getSummary}
\alias{getSummary}
\title{Compute summary results}
\description{
  Compute summary results and the optimal rank from the matrix containing all results                        
}
\usage{ getSummary(obj, NC, NR=96)
}
\arguments{
  \item{obj}{Matrix containing all results in the return list from \code{\link{suitor}}.}
  \item{NC}{The number of columns in \code{data} when \code{\link{suitor}} was called.}
  \item{NR}{The number of rows in \code{data} when \code{\link{suitor}} was called.
            The default is 96.}
}
\details{
  The input matrix \code{obj} must have column 1 as the rank, column 2 as the value
  of  \code{k} in \code{1:k.fold}, column 4 as the training errors, and column 5 as
  the testing errors.
}
\value{ A list containing the objects:
  \itemize{
    \item \code{rank}: The optimal rank
    \item \code{all.results}: Matrix containing training and testing errors for all values
                             of seeds, ranks, folds. NA values appear for runs in which the EM 
                             algorithm did not converge.
    \item{summary}: Data frame of summarized results for each possible rank created from \code{all.results}.
                    The \code{MSErr} column is defined as
                    sqrt(\{fold1 + ... +foldK\}/\{nrow(data)*ncol(data)\})
  }
}
\references{
}
\author{
Donghyuk Lee <donghyuk.lee@nih.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="SUITOR")
  ret <- getSummary(results, 30)
  ret$summary
  ret$rank
}
\keyword{ NMF, cross-validation, mutational signatures } 
\seealso{ \code{\link{plotErrors}} }
