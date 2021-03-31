#' @name operators
#' @title Mathematical operators on CSR matrices
#' @description Implements some mathematical operators between CSR matrices
#' (a.k.a. RsparseMatrix) such as addition and multiplication,
#'  and between CSR matrices and numeric constants, such as division and multiplication.
#' @details The indices of the matrices might be sorted in-place for some operations
#' (see \link{sort_sparse_indices}).
#' @param e1 A CSR matrix.
#' @param e2 Right-hand side of the operation. Either a CSR matrix or a scalar
#' depending on the specific operation. If the RHS is not covered by the functions
#' in this package, will default to the operators from `Matrix` (e.g. when the
#' RHS is a vector)
#' @return A CSR matrix of class `dgRMatrix` (`lgRMatrix` for some of the logical
#' operators).
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- as.csr.matrix(rsparsematrix(4, 3, .5))
#' X + X
#' X * 2
#' X ^ 2
#' ### here the result will be CSC
#' X ^ c(1,2)
NULL

multiply_csr_by_csr <- function(e1, e2) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")

    if (inherits(e1, "nsparseMatrix") && inherits(e2, "nsparseMatrix")) {
        if (is_same_ngRMatrix(e1@p, e2@p, e1@j, e2@j))
            return(e1)
    }
    e1 <- as.csr.matrix(e1)
    e2 <- as.csr.matrix(e2)
    e1 <- sort_sparse_indices(e1)
    e2 <- sort_sparse_indices(e2)
    res <- multiply_csr_elemwise(e1@p, e2@p, e1@j, e2@j, e1@x, e2@x)
    out <- new("dgRMatrix")
    out@Dim <- e1@Dim
    out@Dimnames <- e1@Dimnames
    out@p <- res$indptr
    out@j <- res$indices
    out@x <- res$values
    return(out)
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="sparseMatrix"), multiply_csr_by_csr)

#' @rdname operators
#' @export
setMethod("*", signature(e1="ngRMatrix", e2="sparseMatrix"), multiply_csr_by_csr)

#' @rdname operators
#' @export
setMethod("*", signature(e1="lgRMatrix", e2="sparseMatrix"), multiply_csr_by_csr)

multiply_csr_by_dense <- function(e1, e2) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")

    e1 <- as.csr.matrix(e1)
    if (typeof(e2) == "double") {
        res <- multiply_csr_by_dense_elemwise_double(e1@p, e1@j, e1@x, e2)
    } else if (typeof(e2) == "integer") {
        res <- multiply_csr_by_dense_elemwise_int(e1@p, e1@j, e1@x, e2)
    } else if (typeof(e2) == "logical") {
        res <- multiply_csr_by_dense_elemwise_bool(e1@p, e1@j, e1@x, e2)
    } else {
        mode(e2) <- "double"
        return(e1 * e2)
    }

    out <- e1
    out@x <- res
    return(out)
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="matrix"), multiply_csr_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="ngRMatrix", e2="matrix"), multiply_csr_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="lgRMatrix", e2="matrix"), multiply_csr_by_dense)

add_csr_matrices <- function(e1, e2, is_substraction=FALSE) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to add/substract them.")

    if (inherits(e1, "nsparseMatrix") && inherits(e2, "nsparseMatrix")) {
        if (is_same_ngRMatrix(e1@p, e2@p, e1@j, e2@j)) {
            if (!is_substraction) {
                return(e1)
            } else {
                out <- new("ngRMatrix")
                out@p <- integer(length(e1@p))
                out@Dim <- e1@Dim
                out@Dimnames <- e1@Dimnames
                return(out)
            }
        }
    }
    e1 <- as.csr.matrix(e1)
    e2 <- as.csr.matrix(e2)
    e1 <- sort_sparse_indices(e1)
    e2 <- sort_sparse_indices(e2)
    res <- add_csr_elemwise(e1@p, e2@p, e1@j, e2@j, e1@x, e2@x, is_substraction)
    out <- new("dgRMatrix")
    out@Dim <- e1@Dim
    out@Dimnames <- e1@Dimnames
    out@p <- res$indptr
    out@j <- res$indices
    out@x <- res$values
    return(out)
}

#' @rdname operators
#' @export
setMethod("+", signature(e1="RsparseMatrix", e2="sparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, FALSE))
})

#' @rdname operators
#' @export
setMethod("+", signature(e1="ngRMatrix", e2="sparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, FALSE))
})

#' @rdname operators
#' @export
setMethod("+", signature(e1="lgRMatrix", e2="sparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, FALSE))
})


#' @rdname operators
#' @export
setMethod("-", signature(e1="RsparseMatrix", e2="sparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, TRUE))
})

#' @rdname operators
#' @export
setMethod("-", signature(e1="ngRMatrix", e2="sparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, TRUE))
})

#' @rdname operators
#' @export
setMethod("-", signature(e1="lgRMatrix", e2="sparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, TRUE))
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2)) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 * e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x * as.numeric(e2)
    return(e1)
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2)) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 * e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x * e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2)) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 * e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x * e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("/", signature(e1="RsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 / e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x / as.numeric(e2)
    return(e1)
})

#' @rdname operators
#' @export
setMethod("/", signature(e1="RsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 / e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x / e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("/", signature(e1="RsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 / e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x / e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="RsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %/% e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x %/% as.numeric(e2)
    return(e1)
})

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="RsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %/% e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x %/% e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="RsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %/% e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x %/% e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("%%", signature(e1="RsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %% e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x %% as.numeric(e2)
    return(e1)
})

#' @rdname operators
#' @export
setMethod("%%", signature(e1="RsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %% e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x %% e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("%%", signature(e1="RsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %% e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x %% e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("^", signature(e1="RsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 <= 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 ^ e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x ^ as.numeric(e2)
    return(e1)
})

#' @rdname operators
#' @export
setMethod("^", signature(e1="RsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 <= 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 ^ e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x ^ e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("^", signature(e1="RsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 <= 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 ^ e2)
    }
    e1 <- as.csr.matrix(e1)
    e1@x <- e1@x ^ e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("&", signature(e1="RsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || !isTRUE(e2)) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 & e2)
    }
    e1 <- as.csr.matrix(e1, logical=TRUE)
    e1@x <- e1@x & e2
    return(e1)
})

#' @rdname operators
#' @export
setMethod("|", signature(e1="RsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || typeof(e2) != "logical" || !isFALSE(e2)) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 | e2)
    }
    e1 <- as.csr.matrix(e1, logical=TRUE)
    e1@x <- e1@x | e2
    return(e1)
})
