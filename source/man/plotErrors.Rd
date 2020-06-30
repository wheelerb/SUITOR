
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
              \code{\link{suitor}} or from \code{\link{getSummary}}.}
}
\details{
  The optimal rank is the minimum at which the test error is attained, and appears as
  a red dot on the graph.
}
\value{ NULL
}
\references{
}
\author{
Donghyuk Lee <donghyuk.lee@nih.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="SUITOR")
  s <- getSummary(results, 30)
  plotErrors(s$summary)
}
\keyword{ NMF, cross-validation, mutational signatures } 
%\seealso{ \code{\link{getSummary}} }
