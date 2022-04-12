
\name{SUITOR-package}
\alias{SUITOR}
\docType{package}
\title{
    Number of mutational signatures
}
\description{
    To select the number of mutational signatures through cross-validation.
}
\details{
SUITOR (Selecting the nUmber of mutatIonal signaTures thrOugh 
cRoss-validation), 
an unsupervised cross-validation method that requires little assumptions and 
no numerical approximations to select the optimal number of signatures 
without overfitting the data. 
The full dataset of mutation counts is split into a training set and a 
validation set; for a given number of signatures, these signatures are 
estimated in the training set and then they are used to predict the 
mutations in the validation set. Multiple candidate numbers of signatures 
are considered; and the number of signatures which predicts most closely 
the mutations in the validation set is selected.

The two main functions in this package are \code{\link{suitor}} and 
\code{\link{suitorExtractWH}}.
}
\references{
Lee, D., Wang, D., Yang, X., Shi, J., Landi, M., Zhu, B. (2021) 
SUITOR: selecting the number of mutational signatures through cross-validation.
bioRxiv, doi: https://doi.org/10.1101/2021.07.28.454269. 
}
\author{Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>}
\keyword{package}