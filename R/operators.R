#' @name operators
#' @title Mathematical operators on CSR and COO matrices
#' @description Implements some mathematical operators between CSR
#' (a.k.a. RsparseMatrix) / COO (a.k.a. TsparseMatrix) matrices such as addition and multiplication,
#'  and between CSR/COO matrices and numeric constants, such as division and multiplication.
#' @details The indices of the matrices might be sorted in-place for some operations
#' (see \link{sort_sparse_indices}).
#' 
#' \bold{Important:}: Multiplying NAs by zero here will be treated differently from base R,
#' as it will assume 0*NA = 0 (no entry in the matrix) vs. base R's 0*NA=NA. In order to get the
#' same behavior as in base R, the operations should be done in CSC format.
#' @param e1 A sparse matrix in CSR or COO format, or a scalar or vector, depeding on the operation.
#' @param e2 Another sparse matrix in CSR or COO format, or a scalar or vector, depeding on the operation.
#' @return A CSR or COO matrix depending on the input type and operation.
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- as.csr.matrix(rsparsematrix(4, 3, .5))
#' X + X
#' X * 2
#' X ^ 2
NULL

multiply_csr_by_csr <- function(e1, e2) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")

    if (inherits(e1, "ngRMatrix") && inherits(e2, "ngRMatrix")) {
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

multiply_csr_by_coo <- function(e1, e2) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        warning("Matrices to multiply have different dimensions.")

    e1 <- as.csr.matrix(e1)
    e2 <- as.coo.matrix(e2)
    e1 <- sort_sparse_indices(e1)
    res <- multiply_csr_by_coo_elemwise(
        e1@p, e1@j, e1@x,
        e2@i, e2@j, e2@x,
        nrow(e1), ncol(e1)
    )
    out <- new("dgTMatrix")
    out@i <- res$row
    out@j <- res$col
    out@x <- res$val
    out@Dim <- as.integer(c(max(nrow(e1), nrow(e2)), max(ncol(e1), ncol(e2))))
    return(out)
}

multiply_csr_by_csrorcoo <- function(e1, e2) {
    if (inherits(e1, "RsparseMatrix") && inherits(e2, "RsparseMatrix")) {
        return(multiply_csr_by_csr(e1, e2))
    } else  if (inherits(e1, "RsparseMatrix")) {
        if (inherits(e2, "TsparseMatrix"))
            return(multiply_csr_by_coo(e1, e2))
        else
            return(e1 * as.csr.matrix(e2))
    } else if (inherits(e2, "RsparseMatrix")) {
        if (inherits(e1, "TsparseMatrix"))
            return(multiply_csr_by_coo(e2, e1))
        else
            return(as.csr.matrix(e1) * e2)
    } else {
        e1 <- as.csr.matrix(e1)
        return(multiply_csr_by_csrorcoo(e1, e2))
    }
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="sparseMatrix"), multiply_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("*", signature(e1="ngRMatrix", e2="sparseMatrix"), multiply_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("*", signature(e1="lgRMatrix", e2="sparseMatrix"), multiply_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("*", signature(e1="sparseMatrix", e2="RsparseMatrix"), multiply_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("*", signature(e1="sparseMatrix", e2="ngRMatrix"), multiply_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("*", signature(e1="sparseMatrix", e2="lgRMatrix"), multiply_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("*", signature(e1="CsparseMatrix", e2="TsparseMatrix"), function(e1, e2) {
    return(t(t_shallow(e1) * t(e2)))
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="TsparseMatrix", e2="CsparseMatrix"), function(e1, e2) {
    return(t(t(e1) * t_shallow(e2)))
})

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

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="RsparseMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="ngRMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="lgRMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

multiply_coo_by_dense <- function(e1, e2) {
    
    if (inherits(e2, "float32")) {
        if (is.vector(e2@Data))
            e2@Data <- matrix(e2@Data, ncol=1)
        if (ncol(e2@Data) < ncol(e1)) {
            full_repeats <- ncol(e1) %/% ncol(e2@Data)
            remainder <- ncol(e1) %% ncol(e2@Data)
            e2@Data <- rep(as.vector(e2@Data), full_repeats)
            if (remainder)
                e2@Data <- c(e2@Data, e2@Data[seq(1L, remainder)])
            e2@Data <- matrix(e2@Data, ncol=1)
        }
    }

    if (nrow(e2) < nrow(e1) || ncol(e2) < ncol(e1))
        stop("Cannot multiply matrices - dimensions do not match.")

    if (inherits(e1, "symmetricMatrix") ||
        !inherits(e1, "dsparseMatrix") ||
        (.hasSlot(e1, "diag") && e1@diag != "N")
    ) { 
        e1 <- as.coo.matrix(e1)
    }

    if (typeof(e2) == "double") {
        res <- multiply_coo_by_dense_numeric(
            e2,
            e1@i,
            e1@j,
            e1@x
        )
    } else if (typeof(e2) == "integer") {
        res <- multiply_coo_by_dense_integer(
            e2,
            e1@i,
            e1@j,
            e1@x
        )
    } else if (typeof(e2) == "logical") {
        res <- multiply_coo_by_dense_logical(
            e2,
            e1@i,
            e1@j,
            e1@x
        )
    } else if (inherits(e1, "float32")) {
        res <- multiply_coo_by_dense_float32(
            e2@Data,
            e1@i,
            e1@j,
            e1@x
        )
    } else {
        stop("Internal error. Please open an issue in GitHub.")
    }

    out <- new("dgTMatrix")
    out@i <- res$row
    out@j <- res$col
    out@x <- res$val
    out@Dim <- as.integer(c(max(nrow(e1), nrow(e2)), max(ncol(e1), ncol(e2))))
    return(out)
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="TsparseMatrix", e2="matrix"), multiply_coo_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="TsparseMatrix", e2="float32"), multiply_coo_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="ngTMatrix", e2="matrix"), multiply_coo_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="lgTMatrix", e2="matrix"), multiply_coo_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="TsparseMatrix"), function(e1, e2) multiply_coo_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="TsparseMatrix"), function(e1, e2) multiply_coo_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="ngTMatrix"), function(e1, e2) multiply_coo_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="lgTMatrix"), function(e1, e2) multiply_coo_by_dense(e2, e1))


add_csr_matrices <- function(e1, e2, is_substraction=FALSE) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to add/substract them.")

    if (inherits(e1, "ngRMatrix") && inherits(e2, "ngRMatrix")) {
        if (is_same_ngRMatrix(e1@p, e2@p, e1@j, e2@j)) {
            if (!is_substraction) {
                return(e1)
            } else {
                out <- new("dgRMatrix")
                out@p <- e1@p
                out@j <- e1@j
                out@x <- rep(2., length(e1@j))
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
setMethod("+", signature(e1="sparseMatrix", e2="RsparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e2, e1, FALSE))
})

#' @rdname operators
#' @export
setMethod("+", signature(e1="sparseMatrix", e2="ngRMatrix"), function(e1, e2) {
    return(add_csr_matrices(e2, e1, FALSE))
})

#' @rdname operators
#' @export
setMethod("+", signature(e1="sparseMatrix", e2="lgRMatrix"), function(e1, e2) {
    return(add_csr_matrices(e2, e1, FALSE))
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
setMethod("-", signature(e1="sparseMatrix", e2="RsparseMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, TRUE))
})

#' @rdname operators
#' @export
setMethod("-", signature(e1="sparseMatrix", e2="ngRMatrix"), function(e1, e2) {
    return(add_csr_matrices(e1, e2, TRUE))
})

#' @rdname operators
#' @export
setMethod("-", signature(e1="sparseMatrix", e2="lgRMatrix"), function(e1, e2) {
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
setMethod("*", signature(e1="integer", e2="RsparseMatrix"), function(e1, e2) {
    return(e2 * e1)
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="numeric", e2="RsparseMatrix"), function(e1, e2) {
    return(e2 * e1)
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="logical", e2="RsparseMatrix"), function(e1, e2) {
    return(e2 * e1)
})

### TODO: implement the same ones above but for COO

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
setMethod("/", signature(e1="TsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 / e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("/", signature(e1="TsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 / e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("/", signature(e1="TsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 / e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("%/%", signature(e1="TsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %/% e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("%/%", signature(e1="TsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %/% e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("%/%", signature(e1="TsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %/% e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("%%", signature(e1="TsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %% e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("%%", signature(e1="TsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %% e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("%%", signature(e1="TsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 == 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 %% e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("^", signature(e1="TsparseMatrix", e2="integer"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 <= 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 ^ e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("^", signature(e1="TsparseMatrix", e2="numeric"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 <= 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 ^ e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("^", signature(e1="TsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || e2 <= 0) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 ^ e2)
    }
    e1 <- as.coo.matrix(e1)
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
setMethod("&", signature(e1="TsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || !isTRUE(e2)) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 & e2)
    }
    e1 <- as.coo.matrix(e1, logical=TRUE)
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

#' @rdname operators
#' @export
setMethod("|", signature(e1="TsparseMatrix", e2="logical"), function(e1, e2) {
    if (NROW(e2) != 1L || is.na(e2) || typeof(e2) != "logical" || !isFALSE(e2)) {
        e1 <- as(e1, "CsparseMatrix")
        return(e1 | e2)
    }
    e1 <- as.coo.matrix(e1, logical=TRUE)
    e1@x <- e1@x | e2
    return(e1)
})
