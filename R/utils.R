#' @useDynLib MatrixExtra, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @import methods

#' @title Sort the indices of a sparse matrix or sparse vector
#' @details Will sort the indices of a sparse matrix or sparse vector.
#'
#' In general, the indices of sparse matrices are always meant to be sorted, and
#' it should be rare to have a matrix with unsorted indices when the matrix is the
#' output of some built-in operation from either `Matrix` or `MatrixExtra`, but when
#' the matrices are constructed manually this function can come in handy.
#'
#' \bold{Important:} the input matrix will be modified in-place.
#' @param X A sparse matrix in CSR, CSC, or COO format; or a sparse vector
#' (from the `Matrix` package.)
#' @return The same input `X` with its indices sorted, as invisible (no auto print).
#' Note that the input is itself modified, so there is no need to reassign it.
#' @export
sort_sparse_indices <- function(X) {
    if (inherits(X, "RsparseMatrix")) {

        if (inherits(X, "dsparseMatrix")) {
            sort_sparse_indices_numeric(
                X@p,
                X@j,
                X@x
            )
        } else if (inherits(X, "lsparseMatrix")) {
            sort_sparse_indices_logical(
                X@p,
                X@j,
                X@x
            )
        } else if (inherits(X, "nsparseMatrix")) {
            sort_sparse_indices_binary(
                X@p,
                X@j
            )
        } else {
            X <- as.csr.matrix(X)
            return(sort_sparse_indices(X))
        }

    } else if (inherits(X, "CsparseMatrix")) {

        if (inherits(X, "dsparseMatrix")) {
            sort_sparse_indices_numeric(
                X@p,
                X@i,
                X@x
            )
        } else if (inherits(X, "lsparseMatrix")) {
            sort_sparse_indices_logical(
                X@p,
                X@i,
                X@x
            )
        } else if (inherits(X, "nsparseMatrix")) {
            sort_sparse_indices_binary(
                X@p,
                X@i
            )
        } else {
            X <- as.csc.matrix(X)
            return(sort_sparse_indices(X))
        }

    } else if (inherits(X, "TsparseMatrix")) {

        if (inherits(X, "dsparseMatrix")) {
            sort_coo_indices_numeric(
                X@i,
                X@j,
                X@x
            )
        } else if (inherits(X, "lsparseMatrix")) {
            sort_coo_indices_logical(
                X@i,
                X@j,
                X@x
            )
        } else if (inherits(X, "nsparseMatrix")) {
            sort_coo_indices_binary(
                X@i,
                X@j
            )
        } else {
            X <- as.coo.matrix(X)
            return(sort_sparse_indices(X))
        }

    } else if (inherits(X, "sparseVector")) {

        if (inherits(X, "dsparseVector")) {
            sort_vector_indices_numeric(
                X@i,
                X@x
            )
        } else if (inherits(X, "isparseVector")) {
            sort_vector_indices_integer(
                X@i,
                X@x
            )
        } else if (inherits(X, "lsparseVector")) {
            sort_vector_indices_logical(
                X@i,
                X@x
            )
        } else if (inherits(X, "nsparseVector")) {
            sort_vector_indices_binary(
                X@i
            )
        } else {
            X <- as(X, "dsparseVector")
            return(sort_sparse_indices(X))
        }

    } else {
        stop("Method is only applicable to sparse matrices in CSR, CSC, and COO formats, and to sparse vectors.")
    }
    return(invisible(X))
}
