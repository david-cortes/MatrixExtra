#' @useDynLib MatrixExtra, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @import methods

#' @title MatrixExtra
#' @description Additional methods for the sparse matrices and sparse vector classes
#' from the `Matrix` package, with an emphasis on the CSR (compressed sparse row) format
#' (a.k.a. "RsparseMatrix" in `Matrix` parlance).
#'
#' This package provides, among others:\itemize{
#' \item Fast and multi-threaded matrix multiplications for various combinations of
#' `<sparse,dense>` and `<dense,sparse>` types, including those from the `float` package
#' (see \link[MatrixExtra]{matmult}).
#' \item Operations that work efficiently in the CSR format, such as concatenating by rows
#' (see \link[MatrixExtra]{rbind2-method}) or slicing (see \link[MatrixExtra]{slice}).
#' \item Convenience conversion functions for the different sparse matrix types
#' (see \link[MatrixExtra]{conversions}).
#' \item Overloaded operators such as "+" for `<RsparseMatrix,RsparseMatrix>`
#' and `<RsparseMatrix,scalar>` (see \link[MatrixExtra]{operators}).
#' \item Overloaded mathematical functions that work only on the non-zero entries of matrices
#' (see \link[MatrixExtra]{mathematical-functions}).
#' \item Fast transposes which change the storage format (see \link[MatrixExtra]{t_shallow}).
#' \item Utility functions for sparse matrices (see e.g.
#' \link[MatrixExtra]{sort_sparse_indices}).
#' }
#' 
#' \bold{Important:} `MatrixExtra` changes some important behaviors from the `Matrix` library.
#' See \link{MatrixExtra-options}.
#' @docType package
#' @name MatrixExtra
#' @keywords package
NULL
