#' @importClassesFrom float float32
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads

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
#' @details Be aware that sparse-dense matrix multiplications might suffer from reduced
#' numerical precision, especially when using objects of type `float32`
#' (from the `float` package).
#'
#' When multiplying a sparse matrix by a sparse vector, their indices
#' will be sorted in-place (see \link{sort_sparse_indices}).
#' @param x,y dense (\code{matrix} / \code{float32})
#' and sparse (\code{RsparseMatrix} / \code{CsparseMatrix}) matrices.
#' @return
#' A dense \code{matrix} object.
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

### TODO: one of these float32-by-csr can be simplified to the matrix-vector product version

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="matrix", y="CsparseMatrix"), function(x, y) {
    check_dimensions_match(x, y, matmult=TRUE)

    # restore on exit
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))

    # set num threads to 1 in order to avoid thread contention between BLAS and openmp threads
    RhpcBLASctl::blas_set_num_threads(1L)

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
    if (is.vector(x@Data))
        x@Data <- matrix(x@Data, ncol=1L)
    check_dimensions_match(x, y, matmult=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()

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
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    RhpcBLASctl::blas_set_num_threads(1L)

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
    if (is.vector(x@Data))
        x@Data <- matrix(x@Data, ncol=1L)
    check_dimensions_match(x, y, tcrossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()

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
    check_dimensions_match(x, y, crossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    RhpcBLASctl::blas_set_num_threads(1L)

    if (typeof(x) != "double") mode(x) <- "double"
    y <- as.csc.matrix(y)

    res <- crossprod_dense_csc_numeric(
        x,
        y@p,
        y@i,
        y@x,
        nthreads
    )
    res <- set_dimnames(res, x, y, crossprod=TRUE)
    return(res)
})

#' @rdname matmult
#' @export
setMethod("crossprod", signature(x="float32", y="CsparseMatrix"), function(x, y) {
    if (is.vector(x@Data))
        x@Data <- matrix(x@Data, ncol=1L)
    check_dimensions_match(x, y, crossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()

    y <- as.csc.matrix(y)

    res <- crossprod_dense_csc_float32(
        x@Data,
        y@p,
        y@i,
        y@x,
        nthreads
    )
    res <- set_dimnames(res, x, y, crossprod=TRUE)
    return(new("float32", Data=res))
})

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="matrix"), function(x, y) {
    check_dimensions_match(x, y, matmult=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    RhpcBLASctl::blas_set_num_threads(1L)

    if (typeof(y) != "double") mode(y) <- "double"
    x <- as.csr.matrix(x)

    res <- matmul_csr_dense_numeric(
        x@p,
        x@j,
        x@x,
        y,
        nthreads
    )

    res <- set_dimnames(res, x, y, matmult=TRUE)
    return(res)
})

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="float32"), function(x, y) {
    if (is.vector(y@Data))
        return(gemv_csr_vec(x, y))

    check_dimensions_match(x, y, matmult=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()

    x <- as.csr.matrix(x)

    res <- matmul_csr_dense_float32(
        x@p,
        x@j,
        x@x,
        y@Data,
        nthreads
    )

    res <- set_dimnames(res, x, y, matmult=TRUE)
    return(new("float32", Data=res))
})

#' @rdname matmult
#' @export
setMethod("tcrossprod", signature(x="RsparseMatrix", y="matrix"), function(x, y) {
    check_dimensions_match(x, y, tcrossprod=TRUE)
    nthreads <- RhpcBLASctl::blas_get_num_procs()
    on.exit(RhpcBLASctl::blas_set_num_threads(nthreads))
    RhpcBLASctl::blas_set_num_threads(1L)

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
        } else if (typeof(y) == "float32") {
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

    if (!is.null(rownames(x))) rownames(res) <- rownames(x)
    if (typeof(y) == 'float32') res <- new("float32", Data=res)
    return(matrix(res, ncol=1))
}

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="numeric"), gemv_csr_vec)

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="logical"), gemv_csr_vec)

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="integer"), gemv_csr_vec)

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="RsparseMatrix", y="sparseVector"), gemv_csr_vec)
