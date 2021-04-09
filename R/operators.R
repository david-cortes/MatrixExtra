#' @name operators
#' @title Mathematical operators on CSR and COO matrices
#' @description Implements some mathematical operators between CSR
#' (a.k.a. RsparseMatrix), COO (a.k.a. TsparseMatrix), and dense matrices, such as addition and multiplication,
#'  and between CSR/COO matrices and numeric constants, such as division and multiplication.
#' @details By default, when doing elementwise multiplication (`*`) between a sparse and a dense
#' matrix or vice-versa, if the dense matrix has missing values (`NA` / `NaN`) at some coordinate in
#' which the sparse matrix has no present entry, the resulting output will not have an entry there either,
#' which differs from the behavior of `Matrix` and base R, but makes the operation much faster.
#' 
#' If such missing values are to be preserved, this behavior can be changed through the package
#' options (i.e. `options("MatrixExtra.ignore_na" = FALSE)` - see \link{MatrixExtra-options}).
#' 
#' The indices of the matrices might be sorted in-place for some operations
#' (see \link{sort_sparse_indices}).
#' @param e1 A sparse matrix in CSR or COO format, or a scalar or vector, depeding on the operation.
#' @param e2 Another sparse matrix in CSR or COO format, or a scalar or vector, depeding on the operation.
#' @return A CSR or COO matrix depending on the input type and operation.
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(4, 3, .5, repr="R")
#' X + X
#' X * 2
#' X ^ 2
NULL

### TODO: careful when it sort the indices about not rendering the input unusable:
### e.g. if the input was 'lgRMatrix', gets converted to 'dgRMatrix', and then sorted,
### it will become unusable later.

### TODO: make another option for whether to add back missing values when multiplying
### a CSR by dense, as currently the missing values in the dense one would be ommited
### if they are sparse in the CSR one.

multiply_csr_by_csr <- function(e1, e2, logical=FALSE) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")

    if (inherits(e1, "ngRMatrix") && inherits(e2, "ngRMatrix")) {
        if (is_same_ngRMatrix(e1@p, e2@p, e1@j, e2@j))
            return(e1)
    }

    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)

    check_valid_matrix(e1)
    if (inplace_sort)
        e1 <- deepcopy_before_sort(e1, logical=logical)
    e1 <- as.csr.matrix(e1, logical=logical)
    e1 <- sort_sparse_indices(e1, copy=!inplace_sort)

    check_valid_matrix(e2)
    if (inplace_sort)
        e2 <- deepcopy_before_sort(e2, logical=logical)
    e2 <- as.csr.matrix(e2, logical=logical)
    e2 <- sort_sparse_indices(e2, copy=!inplace_sort)

    if (!logical) {
        res <- multiply_csr_elemwise(e1@p, e2@p, e1@j, e2@j, e1@x, e2@x)
        out <- new("dgRMatrix")
    } else {
        res <- logicaland_csr_elemwise(e1@p, e2@p, e1@j, e2@j, e1@x, e2@x)
        out <- new("lgRMatrix")
    }
    out@Dim <- e1@Dim
    out@Dimnames <- e1@Dimnames
    out@p <- res$indptr
    out@j <- res$indices
    out@x <- res$values
    return(out)
}

multiply_csr_by_coo <- function(e1, e2, logical=FALSE) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        warning("Matrices to multiply have different dimensions.")

    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)

    check_valid_matrix(e1)
    if (inplace_sort)
        e1 <- deepcopy_before_sort(e1, logical=logical)
    e1 <- as.csr.matrix(e1, logical=logical)
    e1 <- sort_sparse_indices(e1, copy=!inplace_sort)
    
    e2 <- as.coo.matrix(e2, logical=logical)
    check_valid_matrix(e2)
    
    if (!logical) {
        res <- multiply_csr_by_coo_elemwise(
            e1@p, e1@j, e1@x,
            e2@i, e2@j, e2@x,
            nrow(e1), ncol(e1)
        )
        out <- new("dgTMatrix")
    } else {
        res <- logicaland_csr_by_coo_elemwise(
            e1@p, e1@j, e1@x,
            e2@i, e2@j, e2@x,
            nrow(e1), ncol(e1)
        )
        out <- new("lgTMatrix")
    }
    out@i <- res$row
    out@j <- res$col
    out@x <- res$val
    out@Dim <- as.integer(c(max(nrow(e1), nrow(e2)), max(ncol(e1), ncol(e2))))
    return(out)
}

multiply_csr_by_csrorcoo_internal <- function(e1, e2, logical=FALSE) {
    if (inherits(e1, "RsparseMatrix") && inherits(e2, "RsparseMatrix")) {
        return(multiply_csr_by_csr(e1, e2, logical=logical))
    } else  if (inherits(e1, "RsparseMatrix")) {
        if (inherits(e2, "TsparseMatrix"))
            return(multiply_csr_by_coo(e1, e2, logical=logical))
        else
            return(multiply_csr_by_csr(e1, as.csr.matrix(e2, logical=logical), logical=logical))
    } else if (inherits(e2, "RsparseMatrix")) {
        if (inherits(e1, "TsparseMatrix"))
            return(multiply_csr_by_coo(e2, e1, logical=logical))
        else
            return(multiply_csr_by_csr(as.csr.matrix(e1, logical=logical), e2, logical=logical))
    } else {
        e1 <- as.csr.matrix(e1)
        return(multiply_csr_by_csrorcoo_internal(e1, e2, logical=logical))
    }
}

multiply_csr_by_csrorcoo <- function(e1, e2) {
    return(multiply_csr_by_csrorcoo_internal(e1, e2, FALSE))
}

logicaland_csr_by_csrorcoo <- function(e1, e2) {
    return(multiply_csr_by_csrorcoo_internal(e1, e2, TRUE))
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
    return(t_shallow(multiply_csr_by_csrorcoo(t_shallow(e1), t_shallow(e2))))
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="TsparseMatrix", e2="CsparseMatrix"), function(e1, e2) {
    return(t_shallow(multiply_csr_by_csrorcoo(t_shallow(e1), t_shallow(e2))))
})

#' @rdname operators
#' @export
setMethod("&", signature(e1="RsparseMatrix", e2="sparseMatrix"), logicaland_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("&", signature(e1="ngRMatrix", e2="sparseMatrix"), logicaland_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("&", signature(e1="lgRMatrix", e2="sparseMatrix"), logicaland_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("&", signature(e1="sparseMatrix", e2="RsparseMatrix"), logicaland_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("&", signature(e1="sparseMatrix", e2="ngRMatrix"), logicaland_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("&", signature(e1="sparseMatrix", e2="lgRMatrix"), logicaland_csr_by_csrorcoo)

#' @rdname operators
#' @export
setMethod("&", signature(e1="CsparseMatrix", e2="TsparseMatrix"), function(e1, e2) {
    return(t_shallow(logicaland_csr_by_csrorcoo(t_shallow(e1), t_shallow(e2))))
})

#' @rdname operators
#' @export
setMethod("&", signature(e1="TsparseMatrix", e2="CsparseMatrix"), function(e1, e2) {
    return(t_shallow(logicaland_csr_by_csrorcoo(t_shallow(e1), t_shallow(e2))))
})

recycle_float32_vector <- function(e1, e2) {
    if (inherits(e2, "float32") && is.vector(e2@Data)) {
        e2@Data <- matrix(e2@Data, ncol=1L)
        if (ncol(e2@Data) < ncol(e1)) {
            full_repeats <- ncol(e1) %/% ncol(e2@Data)
            remainder <- ncol(e1) %% ncol(e2@Data)
            e2@Data <- rep(as.vector(e2@Data), full_repeats)
            if (remainder)
                e2@Data <- c(e2@Data, e2@Data[seq(1L, remainder)])
            e2@Data <- matrix(e2@Data, ncol=1)
        }
    }
    return(e2)
}

multiply_csr_by_dense_internal <- function(e1, e2, logical=FALSE) {
    e2 <- recycle_float32_vector(e1, e2)
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")

    keep_NAs <- !logical && !getOption("MatrixExtra.ignore_na", default=FALSE)
    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)

    check_valid_matrix(e1)
    if (keep_NAs && inplace_sort)
        e1 <- deepcopy_before_sort(e1, logical=logical)
    e1 <- as.csr.matrix(e1, logical=logical)
    if (keep_NAs)
        e1 <- sort_sparse_indices(e1, copy=!inplace_sort)

    if (!logical) {
        if (typeof(e2) == "double") {
            res <- multiply_csr_by_dense_elemwise_double(e1@p, e1@j, e1@x, e2)
        } else if (typeof(e2) == "integer") {
            res <- multiply_csr_by_dense_elemwise_int(e1@p, e1@j, e1@x, e2)
        } else if (typeof(e2) == "logical") {
            res <- multiply_csr_by_dense_elemwise_bool(e1@p, e1@j, e1@x, e2)
        } else if (inherits(e2, "float32")) {
            res <- multiply_csr_by_dense_elemwise_float32(e1@p, e1@j, e1@x, e2@Data)
        } else {
            mode(e2) <- "double"
            return(multiply_csr_by_dense_internal(e1, e2))
        }
    }

    else {
        mode(e2) <- "logical"
        res <- logicaland_csr_by_dense_cpp(e1@p, e1@j, e1@x, e2)
    }

    out <- e1
    out@x <- res

    ### TODO: can this be done more efficiently?
    if (keep_NAs) {

        if (!logical) {

            if (typeof(e2) == "double") {
                res <- add_NAs_from_dense_after_elemenwise_mult_numeric(e1@p, e1@j, e2)
            } else if (typeof(e2) == "integer") {
                res <- add_NAs_from_dense_after_elemenwise_mult_integer(e1@p, e1@j, e2)
            } else if (inherits(e2, "float32")) {
                res <- add_NAs_from_dense_after_elemenwise_mult_float32(e1@p, e1@j, e2@Data)
            }
            
        } else {
            res <- add_NAs_from_dense_after_elemenwise_mult_logical(e1@p, e1@j, e2)
        }

        if (length(res)) {
            out <- as.coo.matrix(out, logical=inherits(out, "lsparseMatrix"))
            X_attr <- attributes(out)
            X_attr$i <- c(X_attr$i, res$ii)
            X_attr$j <- c(X_attr$j, res$jj)
            X_attr$x <- c(X_attr$x, res$xx)
            attributes(out) <- X_attr
        }
    }
    return(out)
}

multiply_csr_by_dense <- function(e1, e2) {
    return(multiply_csr_by_dense_internal(e1, e2, FALSE))
}

logicaland_csr_by_dense <- function(e1, e2) {
    return(multiply_csr_by_dense_internal(e1, e2, TRUE))
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
setMethod("*", signature(e1="RsparseMatrix", e2="float32"), multiply_csr_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="ngRMatrix", e2="float32"), multiply_csr_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="lgRMatrix", e2="float32"), multiply_csr_by_dense)


#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="RsparseMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="ngRMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="lgRMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="RsparseMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="ngRMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="lgRMatrix"), function(e1, e2) multiply_csr_by_dense(e2, e1))


#' @rdname operators
#' @export
setMethod("&", signature(e1="RsparseMatrix", e2="matrix"), logicaland_csr_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="ngRMatrix", e2="matrix"), logicaland_csr_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="lgRMatrix", e2="matrix"), logicaland_csr_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="matrix", e2="RsparseMatrix"), function(e1, e2) logicaland_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("&", signature(e1="matrix", e2="ngRMatrix"), function(e1, e2) logicaland_csr_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("&", signature(e1="matrix", e2="lgRMatrix"), function(e1, e2) logicaland_csr_by_dense(e2, e1))


multiply_coo_by_dense_internal <- function(e1, e2, logical=FALSE) {

    if (!logical && !getOption("MatrixExtra.ignore_na", default=FALSE) && anyNA(e2)) {
        e1 <- as.csc.matrix(e1, logical=inherits(e1, "lsparseMatrix"), binary=inherits(e1, "nsparseMatrix"))
        return(multiply_csc_by_dense_internal(e1, e2, logical))
    }
    
    e2 <- recycle_float32_vector(e1, e2)
    if (nrow(e2) < nrow(e1) || ncol(e2) < ncol(e1))
        stop("Cannot multiply matrices elementwise - dimensions do not match.")

    if (inherits(e1, "symmetricMatrix") ||
        !inherits(e1, ifelse(logical, "lsparseMatrix", "dsparseMatrix")) ||
        (.hasSlot(e1, "diag") && e1@diag != "N")
    ) { 
        e1 <- as.coo.matrix(e1, logical=logical)
    }

    check_valid_matrix(e1)

    if (!logical) {

        out <- new("dgTMatrix")

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
            throw_internal_error()
        }

    } else {

        out <- new("lgTMatrix")
        mode(e2) <- "logical"
        res <- logicaland_coo_by_dense_logical(
            e2,
            e1@i,
            e1@j,
            e1@x
        )

    }
    out@i <- res$row
    out@j <- res$col
    out@x <- res$val
    out@Dim <- as.integer(c(max(nrow(e1), nrow(e2)), max(ncol(e1), ncol(e2))))
    return(out)
}

multiply_coo_by_dense <- function(e1, e2) {
    return(multiply_coo_by_dense_internal(e1, e2, FALSE))
}

logicaland_coo_by_dense <- function(e1, e2) {
    return(multiply_coo_by_dense_internal(e1, e2, TRUE))
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
setMethod("*", signature(e1="ngTMatrix", e2="float32"), multiply_coo_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="lgTMatrix", e2="float32"), multiply_coo_by_dense)

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

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="ngTMatrix"), function(e1, e2) multiply_coo_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="lgTMatrix"), function(e1, e2) multiply_coo_by_dense(e2, e1))


#' @rdname operators
#' @export
setMethod("&", signature(e1="TsparseMatrix", e2="matrix"), logicaland_coo_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="ngTMatrix", e2="matrix"), logicaland_coo_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="lgTMatrix", e2="matrix"), logicaland_coo_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="matrix", e2="TsparseMatrix"), function(e1, e2) logicaland_coo_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("&", signature(e1="matrix", e2="ngTMatrix"), function(e1, e2) logicaland_coo_by_dense(e2, e1))

#' @rdname operators
#' @export
setMethod("&", signature(e1="matrix", e2="lgTMatrix"), function(e1, e2) logicaland_coo_by_dense(e2, e1))

### TODO: this will not handle names correctly, need to add an extra step at the end

multiply_csc_by_dense_internal <- function(e1, e2, logical=FALSE) {
    e2 <- recycle_float32_vector(e1, e2)
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")
    check_valid_matrix(e1)

    ignore_NAs <- logical || getOption("MatrixExtra.ignore_na", default=FALSE)

    if (ignore_NAs) {

        e1 <- as.csc.matrix(e1, logical=logical)
        if (!logical) {
            if (typeof(e2) == "double") {
                res <- multiply_csc_by_dense_ignore_NAs_numeric(e1@p, e1@i, e1@x, e2)
            } else if (typeof(e2) == "integer") {
                res <- multiply_csc_by_dense_ignore_NAs_integer(e1@p, e1@i, e1@x, e2)
            } else if (typeof(e2) == "logical") {
                res <- multiply_csc_by_dense_ignore_NAs_logical(e1@p, e1@i, e1@x, e2)
            } else if (inherits(e2, "float32")) {
                res <- multiply_csc_by_dense_ignore_NAs_float32(e1@p, e1@i, e1@x, e2@Data)
            } else {
                mode(e2) <- "double"
                return(multiply_csc_by_dense_internal(e1, e2))
            }
        }

        else {
            mode(e2) <- "logical"
            res <- logicaland_csc_by_dense_ignore_NAs(e1@p, e1@i, e1@x, e2)
        }

        X_attr <- attributes(e1)
        X_attr$x <- res
        X_attr$i <- deepcopy_int(X_attr$i)
        attributes(e1) <- X_attr
        return(e1)

    } else {

        inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)
        if (inplace_sort)
            e1 <- deepcopy_before_sort(e1, logical=logical)
        e1 <- as.csc.matrix(e1, logical=logical)
        e1 <- sort_sparse_indices(e1, copy=!inplace_sort)

        if (!logical) {

            if (typeof(e2) == "double") {
                res <- multiply_csc_by_dense_keep_NAs_numeric(e1@p, e1@i, e1@x, e2)
            } else if (typeof(e2) == "integer") {
                res <- multiply_csc_by_dense_keep_NAs_integer(e1@p, e1@i, e1@x, e2)
            } else if (typeof(e2) == "logical") {
                res <- multiply_csc_by_dense_keep_NAs_logical(e1@p, e1@i, e1@x, e2)
            } else if (inherits(e2, "float32")) {
                res <- multiply_csc_by_dense_keep_NAs_float32(e1@p, e1@i, e1@x, e2@Data)
            } else {
                mode(e2) <- "double"
                return(multiply_csc_by_dense_internal(e1, e2))
            }

        } else {

            mode(e2) <- "logical"
            res <- logicaland_csc_by_dense_keep_NAs(e1@p, e1@i, e1@x, e2)

        }

        X_attr <- attributes(e1)
        X_attr$p <- res$indptr
        X_attr$i <- res$indices
        X_attr$x <- res$values
        attributes(e1) <- X_attr
        return(e1)
    }
}

multiply_csc_by_dense <- function(e1, e2) {
    return(multiply_csc_by_dense_internal(e1, e2, FALSE))
}

logicaland_csc_by_dense <- function(e1, e2) {
    return(multiply_csc_by_dense_internal(e1, e2, TRUE))
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="CsparseMatrix", e2="matrix"), multiply_csc_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="CsparseMatrix", e2="float32"), multiply_csc_by_dense)

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="CsparseMatrix"), function(e1, e2) {
    return(multiply_csc_by_dense(e2, e1))
})

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="CsparseMatrix"), function(e1, e2) {
    return(multiply_csc_by_dense(e2, e1))
})

#' @rdname operators
#' @export
setMethod("&", signature(e1="CsparseMatrix", e2="matrix"), logicaland_csc_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="CsparseMatrix", e2="float32"), logicaland_csc_by_dense)

#' @rdname operators
#' @export
setMethod("&", signature(e1="matrix", e2="CsparseMatrix"), function(e1, e2) {
    return(logicaland_csc_by_dense(e2, e1))
})

#' @rdname operators
#' @export
setMethod("&", signature(e1="float32", e2="CsparseMatrix"), function(e1, e2) {
    return(logicaland_csc_by_dense(e2, e1))
})

### TODO: add tests for XOR

add_csr_matrices_internal <- function(e1, e2, is_substraction=FALSE,
                                      is_ampersand=FALSE, is_xor=FALSE) {
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to add/substract them.")

    logical <- is_ampersand || is_xor

    if (inherits(e1, "ngRMatrix") && inherits(e2, "ngRMatrix")) {
        if (is_same_ngRMatrix(e1@p, e2@p, e1@j, e2@j)) {
            if (!is_substraction && !is_xor) {
                return(e1)
            } else if (is_xor) {
                out <- new("lgRMatrix")
                out@p <- rep(0L, nrow(e1)+1L)
                out@Dim <- e1@Dim
                out@Dimnames <- e1@Dimnames
                return(out)
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

    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)

    check_valid_matrix(e1)
    if (inplace_sort)
        e1 <- deepcopy_before_sort(e1, logical=logical)
    e1 <- as.csr.matrix(e1, logical=logical)
    e1 <- sort_sparse_indices(e1, copy=!inplace_sort)

    check_valid_matrix(e2)
    if (inplace_sort)
        e2 <- deepcopy_before_sort(e2, logical=logical)
    e2 <- as.csr.matrix(e2, logical=logical)
    e2 <- sort_sparse_indices(e2, copy=!inplace_sort)

    if (!logical) {
        res <- add_csr_elemwise(e1@p, e2@p, e1@j, e2@j, e1@x, e2@x, is_substraction)
        out <- new("dgRMatrix")
    } else if (is_ampersand) {
        res <- logicalor_csr_elemwise(e1@p, e2@p, e1@j, e2@j, e1@x, e2@x, FALSE)
        out <- new("lgRMatrix")
    } else if (is_xor) {
        res <- logicalor_csr_elemwise(e1@p, e2@p, e1@j, e2@j, e1@x, e2@x, TRUE)
        out <- new("lgRMatrix")
    } else {
        throw_internal_error()
    }


    out@Dim <- e1@Dim
    out@Dimnames <- e1@Dimnames
    out@p <- res$indptr
    out@j <- res$indices
    out@x <- res$values
    return(out)
}

add_csr_matrices <- function(e1, e2, is_substraction=FALSE) {
    return(add_csr_matrices_internal(e1, e2, is_substraction, FALSE, FALSE))
}

logicalor_csr_matrices <- function(e1, e2) {
    return(add_csr_matrices_internal(e1, e2, FALSE, TRUE, FALSE))
}

xor_csr_matrices <- function(e1, e2) {
    return(add_csr_matrices_internal(e1, e2, FALSE, FALSE, TRUE))
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
setMethod("+", signature(e1="CsparseMatrix", e2="TsparseMatrix"), function(e1, e2) {
    return(t_shallow(add_csr_matrices(t_shallow(e1), t_shallow(e2), FALSE)))
})

#' @rdname operators
#' @export
setMethod("+", signature(e1="TsparseMatrix", e2="CsparseMatrix"), function(e1, e2) {
    return(t_shallow(add_csr_matrices(t_shallow(e2), t_shallow(e1), FALSE)))
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
setMethod("-", signature(e1="CsparseMatrix", e2="TsparseMatrix"), function(e1, e2) {
    return(t_shallow(add_csr_matrices(t_shallow(e1), t_shallow(e2), TRUE)))
})

#' @rdname operators
#' @export
setMethod("-", signature(e1="TsparseMatrix", e2="CsparseMatrix"), function(e1, e2) {
    return(t_shallow(add_csr_matrices(t_shallow(e2), t_shallow(e1), TRUE)))
})

#' @rdname operators
#' @export
setMethod("|", signature(e1="RsparseMatrix", e2="sparseMatrix"), logicalor_csr_matrices)

#' @rdname operators
#' @export
setMethod("|", signature(e1="ngRMatrix", e2="sparseMatrix"), logicalor_csr_matrices)

#' @rdname operators
#' @export
setMethod("|", signature(e1="lgRMatrix", e2="sparseMatrix"), logicalor_csr_matrices)

#' @rdname operators
#' @export
setMethod("|", signature(e1="sparseMatrix", e2="RsparseMatrix"), logicalor_csr_matrices)

#' @rdname operators
#' @export
setMethod("|", signature(e1="sparseMatrix", e2="ngRMatrix"), logicalor_csr_matrices)

#' @rdname operators
#' @export
setMethod("|", signature(e1="sparseMatrix", e2="lgRMatrix"), logicalor_csr_matrices)

#' @rdname operators
#' @export
setMethod("|", signature(e1="CsparseMatrix", e2="TsparseMatrix"), function(e1, e2) {
    return(t_shallow(logicalor_csr_matrices(t_shallow(e1), t_shallow(e2))))
})

#' @rdname operators
#' @export
setMethod("|", signature(e1="TsparseMatrix", e2="CsparseMatrix"), function(e1, e2) {
    return(t_shallow(logicalor_csr_matrices(t_shallow(e1), t_shallow(e2))))
})

### This didn't work out
# setMethod("xor", signature(e1="RsparseMatrix", e2="sparseMatrix"), xor_csr_matrices)

# setMethod("xor", signature(e1="ngRMatrix", e2="sparseMatrix"), xor_csr_matrices)

# setMethod("xor", signature(e1="lgRMatrix", e2="sparseMatrix"), xor_csr_matrices)

# setMethod("xor", signature(e1="sparseMatrix", e2="RsparseMatrix"), xor_csr_matrices)

# setMethod("xor", signature(e1="sparseMatrix", e2="ngRMatrix"), xor_csr_matrices)

# setMethod("xor", signature(e1="sparseMatrix", e2="lgRMatrix"), xor_csr_matrices)

# setMethod("xor", signature(e1="CsparseMatrix", e2="TsparseMatrix"), function(e1, e2) {
#     return(t_shallow(xor_csr_matrices(t_shallow(e1), t_shallow(e2))))
# })

# setMethod("xor", signature(e1="TsparseMatrix", e2="CsparseMatrix"), function(e1, e2) {
#     return(t_shallow(xor_csr_matrices(t_shallow(e1), t_shallow(e2))))
# })

####

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
