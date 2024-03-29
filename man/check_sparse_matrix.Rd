% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_sparse_matrix}
\alias{check_sparse_matrix}
\title{Check if the underlying data behind a sparse matrix constitutes a valid object}
\usage{
check_sparse_matrix(X, sort = TRUE, remove_zeros = TRUE)
}
\arguments{
\item{X}{A sparse matrix or sparse vector whose underlying arrays are to be checked.}

\item{sort}{Whether to sort the indices of `X` along the way. For this, will make
deep copies of the indices and values so that there are no issues with external
references being updated along the way.}

\item{remove_zeros}{Whether to remove entries in `X` which have a value of zero
but are nevertheless still present in the sparse representation.}
}
\value{
The same matrix or vector `X` (perhaps sorted or with zeros removed depending
on the parameters). If `X` contains data that doesn't make for a valid sparse matrix
(e.g. different number of values and indices), it will throw an error.
}
\description{
Check if the underlying data behind a sparse matrix constitutes a valid object
}
\details{
Makes checks on the data contained in the sparse matrix or sparse vector
object for whether the data constitutes a valid matrix - e.g. indices must not be
negative or larger than the dimensions, index pointer must match with the indices,
etc.

As a short-hand, can also sort the matrix and remove zeros by calling the respective
functions \link{sort_sparse_indices} and \link{remove_sparse_zeros}.

A sparse matrix or sparse vector should never come out with invalid data from
functions from `Matrix` and `MatrixExtra`, with one exception: some `MatrixExtra`
functions might modify indices in-place, which can cause problems if the same
array of indices / values is also used by another matrix or R object elsewhere.

Otherwise, this function is aimed at making checks on matrices that are manually
constructed.
}
