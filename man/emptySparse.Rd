% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{emptySparse}
\alias{emptySparse}
\title{Create Empty Sparse Matrix}
\usage{
emptySparse(nrow = 0L, ncol = 0L, format = "R", dtype = "d")
}
\arguments{
\item{nrow}{Desired number of rows for the matrix.}

\item{ncol}{Desired number of columns for the matrix.}

\item{format}{Storage format for the matrix. Options are:\itemize{
\item "R", which will output a CSR Matrix ("RsparseMatrix").
\item "C", which will output a CSC Matrix ("CsparseMatrix").
\item "T", which will output a COO/triplets Matrix ("TsparseMatrix").
}}

\item{dtype}{Data type for the matrix. Options are:\itemize{
\item "d", which will output a numeric/double type (e.g. "dgRMatrix").
\item "l", which will output a logical/boolean type.
\item "n", which will output a binary type.
}}
}
\value{
A sparse matrix of general type, with the specific class
determined by `format` and `dtype`.
}
\description{
Creates an empty sparse matrix (all values being zeros)
with the requested format and dimensions. This is a faster alternative
to calling `Matrix::Matrix(0, ...)`.
}
\examples{
### This is very fast despite the large dimensions,
### as no data is held in the resulting object
library(MatrixExtra)
X <- emptySparse(nrow=2^20, ncol=2^25, format="T")
}
