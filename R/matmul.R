#' @importClassesFrom float float32
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads

### Peculiarities about R's matrix-by-vector multiplications (as of v4.0.4)
###
### matmul(Mat, vec)
###     -> If 'Mat' has more than one column, 'vec' is a column vector [n,1]
###     -> If 'Mat' has only one column, 'vec' is a row vector [1,n]
### matmul(vec, Mat)
###     -> If 'Mat' has more than one row, 'vec' is a row vector [1,n]
###     -> If 'Mat' has only one row, 'vec' is a column vector [n,1]
### matmul(vec, vec) -> LHS is a row vector [1,n], RHS is a column vector [n,1]
###
### crossprod(Mat, vec)
###     -> If 'Mat' has more than one row, 'vec' is a column vector [n,1]
###     -> If 'Mat' has one column, 'vec' is a row vector [1,n]
### crossprod(vec, Mat) -> 'vec' is a column vector [n,1]
### crossprod(vec, vec) -> 'vec' is a column vector [n,1]
###
### tcrossprod(Mat, vec)
###    -> If 'Mat' has more than one row and more than one column, will fail
###    -> If 'Mat' has only one row, 'vec' is a row vector [1,n]
###    -> If 'Mat' has only one column, 'vec' is a column vector [n,1]
### tcrossprod(vec, Mat)
###    -> If 'Mat' has more than one column, 'vec' is a row vector [1,n]
###    -> If 'Mat' has only one column, 'vec' is a column vector [n,1]
### tcrossprod(vec, vec): 'vec' is a column vector [n,1]


### TODO: try to make the multiplications preserve the names

check_dimensions_match <- function(x, y, matmult=FALSE, crossprod=FALSE, tcrossprod=FALSE) {
    if (matmult) {
        inner_x <- ncol(x)
        inner_y <- nrow(y)
    } else if (crossprod) {
        inner_x <- nrow(x)
        inner_y <- nrow(y)
    } else if (tcrossprod) {
        inner_x <- ncol(x)
        inner_y <- ncol(y)
    } else {
        stop("Internal error. Please open a bug report in GitHub.")
    }

    if (inner_x != inner_y)
        stop("Matrix dimensions do not match.")
}

set_dimnames <- function(res, x, y, matmult=FALSE, crossprod=FALSE, tcrossprod=FALSE) {
    if (matmult) {
        rnames <- rownames(x)
        cnames <- colnames(y)
    } else if (crossprod) {
        rnames <- colnames(x)
        cnames <- colnames(y)
    } else if (tcrossprod) {
        rnames <- rownames(x)
        cnames <- rownames(y)
    } else {
        stop("Internal error. Please open a bug report in GitHub.")
    }

    if (!is.null(rnames))
        rownames(res) <- rnames
    if (!is.null(cnames))
        colnames(res) <- cnames
    return(res)
}

#' @title Multithreaded Sparse-Dense Matrix and Vector Multiplications
#' @description Multithreaded <matrix, matrix> multiplications
#' (`\%*\%`, `crossprod`, and `tcrossprod`)
#' and <matrix, vector> multiplications (`\%*\%`),
#' for <sparse, dense> matrix combinations and <sparse, vector> combinations
#' (See signatures for supported combinations).
#'
#' Objects from the `float` package are also supported for some combinations.
#' @details Will try to use the same number of threads that are configured for BLAS.
#' These can be set through the `RhpcBLASctl` package (see
#' \link[RhpcBLASctl]{blas_set_num_threads} and \link[RhpcBLASctl]{blas_get_num_procs}).
#'
#' Be aware that sparse-dense matrix multiplications might suffer from reduced
#' numerical precision, especially when using objects of type `float32`
#' (from the `float` package).
#'
#' Internally, these functions use BLAS level-1 routines, so their speed might depend on
#' the BLAS backend being used (e.g. MKL, OpenBLAS) - that means: they might be quite slow
#' on a default install of R for Windows.
#'
#' When multiplying a sparse matrix by a sparse vector, their indices
#' will be sorted in-place (see \link{sort_sparse_indices}).
#'
#' In order to match exactly with base R's behaviors, when passing vectors to these
#' operators, will assume their shape as follows:\itemize{
#' \item MatMult(Matrix, vector): column vector if the matrix has more than one column
#' or is empty, row vector if the matrix has only one column.
#' \item MatMult(vector, Matrix): row vector if the matrix has more than one row,
#' column vector if the matrix has only one row
#' \item crossprod(Matrix, vector): column vector if the matrix has more than one row,
#' row vector if the matrix has only one row.
#' \item crossprod(vector, Matrix): column vector.
#' \item tcrossprod(Matrix, vector): row vector if the matrix has only one row,
#' column vector if the matrix has only one column, and will throw an error otherwise.
#' \item tcrossprod(vector, Matrix): row vector if the matrix has more than one column,
#' column vector if the matrix has only one column.
#' }
#' @param x,y dense (\code{matrix} / \code{float32})
#' and sparse (\code{RsparseMatrix} / \code{CsparseMatrix}) matrices or vectors
#' (\code{sparseVector}, \code{numeric}, \code{integer}, \code{logical}).
#' @return
#' A dense \code{matrix} object in most cases, except for `<RsparseMatrix, sparseVector>`,
#' which will return a CSC matrix (`dgCMatrix`).
#'
#' @name matmult
#' @rdname matmult
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' ## Will use the same number of threads as for BLAS
#' curr_nthreads <- RhpcBLASctl::blas_get_num_procs()
#' ## Set the number of threads here (will restore later)
#' RhpcBLASctl::blas_set_num_threads(1L)
#'
#' ## Generate random matrices
#' set.seed(1)
#' A <- rsparsematrix(5,4,.5)
#' B <- rsparsematrix(4,3,.5)
#'
#' ## Now multiply in some supported combinations
#' as.matrix(A) %*% as.csc.matrix(B)
#' as.csr.matrix(A) %*% as.matrix(B)
#' crossprod(as.matrix(B), as.csc.matrix(B))
#' tcrossprod(as.csr.matrix(A), as.matrix(A))
#'
#' ## Restore the number of threads for BLAS
#' RhpcBLASctl::blas_set_num_threads(curr_nthreads)
NULL

#### Matrices ----

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="matrix", y="CsparseMatrix"), function(x, y) {
    check_dimensions_match(x, y, matmult=TRUE)

    # restore on exit
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    nthreads <- max(nthreads, 1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))

    # set num threads to 1 in order to avoid thread contention between BLAS and openmp threads
    if (nthreads > 1) RhpcBLASctl::blas_set_num_threads(1L)

    if (typeof(x) != "double") mode(x) <- "double"
    y <- as.csc.matrix(y)

    res <- matmul_dense_csc_numeric(
        x,
        y@p,
        y@i,
        y@x,
        nthreads
    )

    res <- set_dimnames(res, x, y, matmult=TRUE)
    return(res)
})

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="float32", y="CsparseMatrix"), function(x, y) {

    if (is.vector(x@Data)) {

        if (inherits(y, "symmetricMatrix") ||
            (.hasSlot(y, "diag") && y@diag != "N") ||
            (.hasSlot(y, "x") && !inherits(y, "dsparseMatrix"))
        ) {
            y <- as.csc.matrix(y)
        }

        ### To match base R, if 'y' has more than one row, 'x' is [1,n], otherwise [n,1]
        if (nrow(y) == 1L) {
            if (.hasSlot(y, "x")) {
                res <- matmul_colvec_by_srowvecascsc(
                    x@Data,
                    y@p,
                    y@i,
                    y@x
                )
            } else {
                res <- matmul_colvec_by_srowvecascsc_binary(
                    x@Data,
                    y@p,
                    y@i
                )
            }
            return(new("float32", Data=res))
        } else {
            if (nrow(y) != length(x@Data))
                stop("(row) vector-Matrix multiplication dimensions do not match.")

            if (.hasSlot(y, "x")) {
                res <- matmul_rowvec_by_csc(
                    x@Data,
                    y@p,
                    y@i,
                    y@x
                )
            } else {
                res <- matmul_rowvec_by_cscbin(
                    x@Data,
                    y@p,
                    y@i
                )
            }
            return(new("float32", Data=res))
        }
    }

    check_dimensions_match(x, y, matmult=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    nthreads <- max(nthreads, 1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    if (nthreads > 1) RhpcBLASctl::blas_set_num_threads(1L)

    y <- as.csc.matrix(y)

    res <- matmul_dense_csc_float32(
        x@Data,
        y@p,
        y@i,
        y@x,
        nthreads
    )

    res <- set_dimnames(res, x, y, matmult=TRUE)
    return(new("float32", Data=res))
})

#' @rdname matmult
#' @export
setMethod("tcrossprod", signature(x="matrix", y="RsparseMatrix"), function(x, y) {
    check_dimensions_match(x, y, tcrossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    nthreads <- max(nthreads, 1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    if (nthreads > 1) RhpcBLASctl::blas_set_num_threads(1L)

    if (typeof(x) != "double") mode(x) <- "double"
    y <- as.csr.matrix(y)

    res <- tcrossprod_dense_csr_numeric(
        x,
        y@p,
        y@j,
        y@x,
        nthreads, ncol(y)
    )
    res <- set_dimnames(res, x, y, tcrossprod=TRUE)
    return(res)
})

#' @rdname matmult
#' @export
setMethod("tcrossprod", signature(x="float32", y="RsparseMatrix"), function(x, y) {

    if (is.vector(x@Data)) {

        if (inherits(y, "symmetricMatrix") ||
            (.hasSlot(y, "diag") && y@diag != "N") ||
            (.hasSlot(y, "x") && !inherits(y, "dsparseMatrix"))
        ) {
            y <- as.csr.matrix(y)
        }

        ### To match with base R, if 'y' has only one column, x is [n,1], otherwise [1,n]
        if (ncol(y) == 1L) {
            if (.hasSlot(y, "x")) {
                res <- matmul_colvec_by_srowvecascsc(
                    x@Data,
                    y@p,
                    y@j,
                    y@x
                )
            } else {
                res <- matmul_colvec_by_srowvecascsc_binary(
                    x@Data,
                    y@p,
                    y@j
                )
            }
            return(new("float32", Data=res))
        } else {
            if (.hasSlot(y, "x")) {
                res <- matmul_rowvec_by_csc(
                    x@Data,
                    y@p,
                    y@j,
                    y@x
                )
            } else {
                res <- matmul_rowvec_by_cscbin(
                    x@Data,
                    y@p,
                    y@j
                )
            }
            return(new("float32", Data=res))
        }
    }

    check_dimensions_match(x, y, tcrossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    nthreads <- max(nthreads, 1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    if (nthreads > 1) RhpcBLASctl::blas_set_num_threads(1L)

    y <- as.csr.matrix(y)

    res <- tcrossprod_dense_csr_float32(
        x@Data,
        y@p,
        y@j,
        y@x,
        nthreads, ncol(y)
    )
    res <- set_dimnames(res, x, y, tcrossprod=TRUE)
    return(new("float32", Data=res))
})


#' @rdname matmult
#' @export
setMethod("crossprod", signature(x="matrix", y="CsparseMatrix"), function(x, y) {
    return(t(x) %*% y)
})

#' @rdname matmult
#' @export
setMethod("crossprod", signature(x="float32", y="CsparseMatrix"), function(x, y) {

    if (is.vector(x@Data)) {
        if (length(x@Data) != nrow(y))
            stop("(column) vector-Matrix crossprod dimensions do not match.")
        if (inherits(y, "symmetricMatrix") ||
            (.hasSlot(y, "diag") && y@diag != "N") ||
            (.hasSlot(y, "x") && !inherits(y, "dsparseMatrix"))
        ) {
            y <- as.csc.matrix(y)
        }
        if (.hasSlot(y, "x")) {
            res <- matmul_rowvec_by_csc(
                x@Data,
                y@p,
                y@i,
                y@x
            )
        } else {
            res <- matmul_rowvec_by_cscbin(
                x@Data,
                y@p,
                y@i
            )
        }
        return(new("float32", Data=res))
    }

    return(t(x) %*% y)
})

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="matrix"), function(x, y) {
    return(tcrossprod(x, t(y)))
})

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="float32"), function(x, y) {

    if (is.vector(y@Data)) {
        if (ncol(x) == 1L) {
            if (!inherits(x, "dsparseMatrix") ||
                inherits(x, "symmetricMatrix") ||
                (.hasSlot(x, "diag") && x@diag != "N")) {
                x <- as.csr.matrix(x)
            }
            res <- matmul_colvec_by_scolvecascsr_f32(
                y@Data,
                x@p,
                x@j,
                x@x
            )
            return(new("float32", Data=res))
        } else {
            return(gemv_csr_vec(x, y))
        }
    }

    return(tcrossprod(x, t(y)))
})

#' @rdname matmult
#' @export
setMethod("tcrossprod", signature(x="RsparseMatrix", y="matrix"), function(x, y) {
    check_dimensions_match(x, y, tcrossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    nthreads <- max(nthreads, 1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    if (nthreads > 1) RhpcBLASctl::blas_set_num_threads(1L)

    if (typeof(y) != "double") mode(y) <- "double"
    x <- as.csr.matrix(x)

    res <- tcrossprod_csr_dense_numeric(
        x@p,
        x@j,
        x@x,
        y,
        nthreads
    )

    res <- set_dimnames(res, x, y, tcrossprod=TRUE)
    return(res)
})

#' @rdname matmult
#' @export
setMethod("tcrossprod", signature(x="RsparseMatrix", y="float32"), function(x, y) {
    check_dimensions_match(x, y, tcrossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    nthreads <- max(nthreads, 1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    if (nthreads > 1) RhpcBLASctl::blas_set_num_threads(1L)

    x <- as.csr.matrix(x)

    res <- tcrossprod_csr_dense_float32(
        x@p,
        x@j,
        x@x,
        y@Data,
        nthreads
    )

    res <- set_dimnames(res, x, y, tcrossprod=TRUE)
    return(new("float32", Data=res))
})

#### Vectors ----

gemv_csr_vec <- function(x, y) {
    if (ncol(x) != length(y))
        stop("Matrix-vector dimensions do not match.")
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    x <- as.csr.matrix(x)

    if (!inherits(y, "sparseVector")) {

        ### dense vectors from base R
        if (typeof(y) == "double") {
            res <- matmul_csr_dvec_numeric(
                x@p,
                x@j,
                x@x,
                y,
                nthreads
            )
        } else if (typeof(y) == "integer") {
            res <- matmul_csr_dvec_integer(
                x@p,
                x@j,
                x@x,
                y,
                nthreads
            )
        } else if (typeof(y) == "logical") {
            res <- matmul_csr_dvec_logical(
                x@p,
                x@j,
                x@x,
                y,
                nthreads
            )
        } else if (inherits(y, "float32")) {
            res <- matmul_csr_dvec_float32(
                x@p,
                x@j,
                x@x,
                y@Data,
                nthreads
            )
        } else {
            y <- as.numeric(y)
            if (typeof(y) != "double")
                mode(y) <- "double"
            return(x %*% y)
        }

    } else {

        ### sparse vectors from matrix
        x <- sort_sparse_indices(x)
        y <- sort_sparse_indices(y)

        if (inherits(y, "dsparseVector")) {
            res <- matmul_csr_svec_numeric(
                x@p,
                x@j,
                x@x,
                y@i,
                y@x,
                nthreads
            )
        } else if (inherits(y, "isparseVector")) {
            res <- matmul_csr_svec_integer(
                x@p,
                x@j,
                x@x,
                y@i,
                y@x,
                nthreads
            )
        } else if (inherits(y, "lsparseVector")) {
            res <- matmul_csr_svec_logical(
                x@p,
                x@j,
                x@x,
                y@i,
                y@x,
                nthreads
            )
        } else if (inherits(y, "nsparseVector")) {
            res <- matmul_csr_svec_binary(
                x@p,
                x@j,
                x@x,
                y@i,
                nthreads
            )
        } else {
            y <- as.numeric(y)
            if (typeof(y) != "double")
                mode(y) <- "double"
            return(x %*% y)
        }
    }

    if (!is.null(rownames(x)))
        names(res) <- rownames(x)

    if (!inherits(y, "float32")) {
        return(matrix(res, ncol=1))
    } else {
        res <- new("float32", Data=matrix(res, ncol=1))
        return(res)
    }
}

outerprod_csrsinglecol_by_dvec <- function(x, y) {
    if (ncol(x) != 1L)
        stop("Internal error. Please create a bug report in GitHub.")

    if (!inherits(x, "dsparseMatrix") ||
        inherits(x, "symmetricMatrix") ||
        (.hasSlot(x, "diag") && x@diag != "N")
    ) {
        x <- as.csr.matrix(x)
    }

    if (inherits(y, "sparseVector")) {
        y <- sort_sparse_indices(y)

        if (inherits(y, "dsparseVector")) {
            res <- matmul_spcolvec_by_scolvecascsr_numeric(
                x@p,
                x@j,
                x@x,
                y@i,
                y@x,
                y@length
            )
        } else if (inherits(y, "isparseVector")) {
            res <- matmul_spcolvec_by_scolvecascsr_integer(
                x@p,
                x@j,
                x@x,
                y@i,
                y@x,
                y@length
            )
        } else if (inherits(y, "lsparseVector")) {
            res <- matmul_spcolvec_by_scolvecascsr_logical(
                x@p,
                x@j,
                x@x,
                y@i,
                y@x,
                y@length
            )
        } else if (inherits(y, "nsparseVector")) {
            res <- matmul_spcolvec_by_scolvecascsr_binary(
                x@p,
                x@j,
                x@x,
                y@i,
                y@length
            )
        } else {
            y <- as(y, "dsparseVector")
            return(outerprod_csrsinglecol_by_dvec(x, y))
        }
        out <- new("dgCMatrix")
        out@p <- res$indptr
        out@i <- res$indices
        out@x <- res$values
        out@Dim <- as.integer(c(nrow(x), y@length))
        out@Dimnames <- list(rownames(x), NULL)
        return(out)
    } else {

        if (typeof(y) != "double")
            mode(y) <- "double"

        res <- matmul_colvec_by_scolvecascsr(
            y,
            x@p,
            x@j,
            x@x
        )
        if (!is.null(rownames(x)))
            rownames(res) <- rownames(x)
        if ("names" %in% names(attributes(y)))
            colnames(x) <- names(y)
        return(res)
    }

}

matmul_csr_vec <- function(x, y) {
    if (ncol(x) == 1L)
        return(outerprod_csrsinglecol_by_dvec(x, y))
    else
        return(gemv_csr_vec(x, y))
}

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="numeric"), matmul_csr_vec)

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="logical"), matmul_csr_vec)

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="integer"), matmul_csr_vec)

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="sparseVector"), matmul_csr_vec)
