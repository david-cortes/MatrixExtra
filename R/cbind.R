#' @name cbind2-method
#' @title Concatenate sparse matrices by columns
#' @description `cbind2` method for the CSR and COO sparse matrix and
#' sparse vector classes from `Matrix`,
#' taking the most efficient route for the concatenation according to the input types.
#' @param x First matrix to concatenate.
#' @param y Second matrix to concatenate.
#' @return A sparse matrix (storage order varying depending on the input types).
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(3, 4, .3)
#' X <- as(X, "TsparseMatrix")
#' inherits(cbind2(X, X), "TsparseMatrix")
NULL

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="TsparseMatrix", y="TsparseMatrix"), function(x, y) {
    return(t_shallow(rbind2_coo(t_shallow(x), t_shallow(y))))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="TsparseMatrix", y="sparseVector"), function(x, y) {
    return(t_shallow(rbind2_coo_vec(t_shallow(x), y, TRUE)))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="sparseVector", y="TsparseMatrix"), function(x, y) {
    return(t_shallow(rbind2_coo_vec(t_shallow(y), x, FALSE)))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="CsparseMatrix", y="sparseVector"), function(x, y) {
    check_valid_matrix(x)
    return(t_shallow(rbind2_generic(t_shallow(x), y)))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="sparseVector", y="CsparseMatrix"), function(x, y) {
    check_valid_matrix(y)
    return(t_shallow(rbind2_generic(x, t_shallow(y))))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="sparseVector", y="sparseVector"), function(x, y) {
    return(t_shallow(rbind2_generic(x, y)))
})


cbind_csr <- function(x, y) {
    check_valid_matrix(x)
    check_valid_matrix(y)

    binary_types <- c("nsparseMatrix", "nsparseVector")
    logical_types <- c("lsparseMatrix", "lsparseVector")
    x_is_binary <- inherits(x, binary_types)
    x_is_logical <- inherits(x, logical_types)
    y_is_binary <- inherits(y, binary_types)
    y_is_logical <- inherits(y, logical_types)

    if (x_is_binary && y_is_binary) {
        if (inherits(x, "symmetricMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
            x <- as.csr.matrix(x, binary=TRUE)
        if (inherits(y, "symmetricMatrix") || (.hasSlot(y, "diag") && y@diag != "N"))
            y <- as.csr.matrix(y, binary=TRUE)
        res <- cbind_csr_binary(x@p, x@j, y@p, y@j + ncol(x))
        out <- new("ngRMatrix")
    } else if ((x_is_binary || x_is_logical) && (y_is_binary || y_is_logical)) {
        x <- as.csr.matrix(x, logical=TRUE)
        y <- as.csr.matrix(y, logical=TRUE)
        res <- cbind_csr_logical(x@p, x@j, x@x, y@p, y@j+ncol(x), y@x)
        out <- new("lgRMatrix")
    } else {
        x <- as.csr.matrix(x)
        y <- as.csr.matrix(y)
        res <- cbind_csr_numeric(x@p, x@j, x@x, y@p, y@j+ncol(x), y@x)
        out <- new("dgRMatrix")
    }

    out@p <- res$indptr
    out@j <- res$indices
    if (.hasSlot(out, "x"))
        out@x <- res$values
    out@Dim <- as.integer(c(max(nrow(x), nrow(y)), ncol(x)+ncol(y)))
    Dimnames <- list(NULL, NULL)
    ### TODO: handle dimnames here
    return(out)
}

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="RsparseMatrix", y="RsparseMatrix"), cbind_csr)

cbind_csr_coo <- function(x, y) {
    if (inherits(x, "TsparseMatrix")) {
        x <- as.csr.matrix(x, logical=inherits(x, "lsparseMatrix"), binary=inherits(x, "nsparseMatrix"))
    } else if (inherits(y, "TsparseMatrix")) {
        y <- as.csr.matrix(y, logical=inherits(y, "lsparseMatrix"), binary=inherits(y, "nsparseMatrix"))
    }
    return(cbind_csr(x, y))
}

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="TsparseMatrix", y="RsparseMatrix"), cbind_csr_coo)

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="RsparseMatrix", y="TsparseMatrix"), cbind_csr_coo)

cbind_csr_vec <- function(x, y) {
    return(t_shallow(rbind2(t_shallow(x), y)))
}

cbind_vec_csr <- function(x, y) {
    return(t_shallow(rbind2(x, t_shallow(y))))
}

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="RsparseMatrix", y="numeric"), cbind_csr_vec)

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="RsparseMatrix", y="integer"), cbind_csr_vec)

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="RsparseMatrix", y="logical"), cbind_csr_vec)

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="RsparseMatrix", y="sparseVector"), cbind_csr_vec)


#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="numeric", y="RsparseMatrix"), cbind_vec_csr)

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="integer", y="RsparseMatrix"), cbind_vec_csr)

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="logical", y="RsparseMatrix"), cbind_vec_csr)

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="sparseVector", y="RsparseMatrix"), cbind_vec_csr)

### TODO: add cbind_csc for batched binding

