% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbind.R
\name{rbind2-method}
\alias{rbind2-method}
\alias{rbind2,RsparseMatrix,RsparseMatrix-method}
\alias{rbind2,sparseVector,RsparseMatrix-method}
\alias{rbind2,RsparseMatrix,sparseVector-method}
\alias{rbind2,sparseVector,sparseVector-method}
\alias{rbind2,CsparseMatrix,CsparseMatrix-method}
\alias{rbind2,sparseVector,CsparseMatrix-method}
\alias{rbind2,CsparseMatrix,sparseVector-method}
\alias{rbind2,RsparseMatrix,CsparseMatrix-method}
\alias{rbind2,CsparseMatrix,RsparseMatrix-method}
\alias{rbind2,RsparseMatrix,TsparseMatrix-method}
\alias{rbind2,TsparseMatrix,RsparseMatrix-method}
\alias{rbind2,RsparseMatrix,numeric-method}
\alias{rbind2,RsparseMatrix,integer-method}
\alias{rbind2,RsparseMatrix,logical-method}
\alias{rbind2,numeric,RsparseMatrix-method}
\alias{rbind2,integer,RsparseMatrix-method}
\alias{rbind2,logical,RsparseMatrix-method}
\alias{rbind2,TsparseMatrix,TsparseMatrix-method}
\alias{rbind2,TsparseMatrix,sparseVector-method}
\alias{rbind2,sparseVector,TsparseMatrix-method}
\alias{rbind2,TsparseMatrix,CsparseMatrix-method}
\alias{rbind2,CsparseMatrix,TsparseMatrix-method}
\title{Concatenate sparse matrices/vectors by rows}
\usage{
\S4method{rbind2}{RsparseMatrix,RsparseMatrix}(x, y)

\S4method{rbind2}{sparseVector,RsparseMatrix}(x, y)

\S4method{rbind2}{RsparseMatrix,sparseVector}(x, y)

\S4method{rbind2}{sparseVector,sparseVector}(x, y)

\S4method{rbind2}{CsparseMatrix,CsparseMatrix}(x, y)

\S4method{rbind2}{sparseVector,CsparseMatrix}(x, y)

\S4method{rbind2}{CsparseMatrix,sparseVector}(x, y)

\S4method{rbind2}{CsparseMatrix,CsparseMatrix}(x, y)

\S4method{rbind2}{RsparseMatrix,CsparseMatrix}(x, y)

\S4method{rbind2}{CsparseMatrix,RsparseMatrix}(x, y)

\S4method{rbind2}{RsparseMatrix,TsparseMatrix}(x, y)

\S4method{rbind2}{TsparseMatrix,RsparseMatrix}(x, y)

\S4method{rbind2}{RsparseMatrix,numeric}(x, y)

\S4method{rbind2}{RsparseMatrix,integer}(x, y)

\S4method{rbind2}{RsparseMatrix,logical}(x, y)

\S4method{rbind2}{numeric,RsparseMatrix}(x, y)

\S4method{rbind2}{integer,RsparseMatrix}(x, y)

\S4method{rbind2}{logical,RsparseMatrix}(x, y)

\S4method{rbind2}{TsparseMatrix,TsparseMatrix}(x, y)

\S4method{rbind2}{TsparseMatrix,sparseVector}(x, y)

\S4method{rbind2}{sparseVector,TsparseMatrix}(x, y)

\S4method{rbind2}{TsparseMatrix,CsparseMatrix}(x, y)

\S4method{rbind2}{CsparseMatrix,TsparseMatrix}(x, y)
}
\arguments{
\item{x}{First matrix to concatenate.}

\item{y}{Second matrix to concatenate.}
}
\value{
A sparse matrix, usually in CSR format but some combinations
might return COO or CSC.
}
\description{
`rbind2` method for the sparse matrix and sparse vector classes from `Matrix`,
taking the most efficient route for the concatenation according to the input types.
}
\examples{
library(Matrix)
library(MatrixExtra)
set.seed(1)
X <- rsparsematrix(3, 4, .3)
X <- as(X, "RsparseMatrix")
inherits(rbind2(X, X), "RsparseMatrix")
inherits(rbind(X, X, as.csc.matrix(X), X), "RsparseMatrix")
inherits(rbind2(as.coo.matrix(X), as.coo.matrix(X)), "TsparseMatrix")
inherits(rbind2(as.csc.matrix(X), as.csc.matrix(X)), "CsparseMatrix")
}
