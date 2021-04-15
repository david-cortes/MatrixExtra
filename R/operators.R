#' @name operators
#' @title Mathematical operators on sparse matrices and sparse vectors
#' @description Implements some mathematical operators between sparse-sparse and sparse-dense
#' matrices and vectors, such as `CSR + CSR`, `CSR + COO`, `CSR * vector`, `CSR * dense`, among
#' others, which typically work natively in the storage order of the inputs without data duplication.
#' @details By default, when doing elementwise multiplication (`*`) between a sparse and a dense
#' matrix or vice-versa, if the dense matrix has missing values (`NA` / `NaN`) at some coordinate in
#' which the sparse matrix has no present entry, the resulting output will not have an entry there either,
#' which differs from the behavior of `Matrix` and base R, but makes the operation much faster.
#' The same applies to division by zero and exponentiation to zero.
#' 
#' If such missing values (or infinites and ones) are to be preserved, this behavior can be
#' changed through the package options (i.e. `options("MatrixExtra.ignore_na" = FALSE)` - see
#' \link{MatrixExtra-options}).
#' 
#' The indices of the matrices might be sorted in-place for some operations
#' (see \link{sort_sparse_indices}).
#' @param e1 A sparse or dense matrix or vector/
#' @param e2 Another sparse or dense matrix or vector.
#' @return A CSR or COO matrix depending on the input type and operation. Some operations
#' (blocked by default) will produce dense matrices as outputs.
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(4, 3, .5, repr="R")
#' X + X
#' X * X
#' X * as.coo.matrix(X)
#' X * 2
#' X * 1:4
#' X ^ 2
#' X ^ (1:4)
#' 
#' ### Beware
#' set_new_matrix_behavior()
#' suppressWarnings(X / 0)
#' restore_old_matrix_behavior()
#' suppressWarnings(X / 0)
NULL

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

### TODO: implement proper recycling of elements
### TODO: the 'allow_remainder' option is no longer used
recycle_float32_vector <- function(e1, e2, allow_remainder=FALSE) {
    if (inherits(e2, "float32") && is.vector(e2@Data)) {
        e2@Data <- matrix(e2@Data, ncol=1L)
        if (ncol(e2@Data) < ncol(e1)) {
            full_repeats <- ncol(e1) %/% ncol(e2@Data)
            remainder <- ncol(e1) %% ncol(e2@Data)
            if (!allow_remainder && remainder > 0)
                stop("Uneven recycling of vector elements not yet supported.")
            e2@Data <- rep(as.vector(e2@Data), full_repeats)
            if (remainder)
                e2@Data <- c(e2@Data, e2@Data[seq(1L, remainder)])
            e2@Data <- matrix(e2@Data, ncol=1)
        }
    }
    return(e2)
}

multiply_csr_by_dense_internal <- function(e1, e2, logical=FALSE) {
    if (inherits(e2, "float32") && is.vector(e2@Data)) {
        if (!NROW(e2@Data)) {
            if (logical)
                return(logical())
            else
                return(numeric())
        }
        if (NROW(e2@Data) > nrow(e1)*ncol(e2))
            stop("Vector to multiply with has more entries than matrix dimensions.")
        return(multiply_csr_by_dense_internal(e1, float::dbl(e2), logical))
    }
    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")

    keep_NAs <- !logical && !getOption("MatrixExtra.ignore_na", default=FALSE)
    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)

    if (keep_NAs || typeof(e2) == "double") {
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, logical=logical, X_is_LHS=TRUE, op="*"))
    }

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

    if (logical || typeof(e2) == "double")
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, logical=logical))

    if (inherits(e2, "float32") && is.vector(e2@Data) && !NROW(e2@Data)) {
        if (logical)
            return(logical())
        else
            return(numeric())
    }

    if (!logical && !getOption("MatrixExtra.ignore_na", default=FALSE) && anyNA(e2)) {
        e1 <- as.csc.matrix(e1, logical=inherits(e1, "lsparseMatrix"), binary=inherits(e1, "nsparseMatrix"))
        return(multiply_csc_by_dense_internal(e1, e2, logical))
    }
    
    e2 <- recycle_float32_vector(e1, e2, FALSE)
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

    ignore_NAs <- logical || getOption("MatrixExtra.ignore_na", default=FALSE)

    if (inherits(e2, "float32") && is.vector(e2@Data)) {
        if (!NROW(e2@Data)) {
            if (logical)
                return(logical())
            else
                return(numeric())
        }
        if (NROW(e2@Data) > nrow(e1)*ncol(e2))
            stop("Vector to multiply with has more entries than matrix dimensions.")
        if (NROW(e2@Data) > nrow(e1) ||
            (NROW(e2@Data) < nrow(e1) && (NROW(e2@Data) %% nrow(e1))) ||
            (NROW(e2@Data) != nrow(e1) && !ignore_NAs && anyNA(e2))
        ) {
            return(e1 * float::dbl(e2))
        } else {
            e2 <- recycle_float32_vector(e1, e2, FALSE)
        }
    }

    if (nrow(e1) != nrow(e2) || ncol(e1) != ncol(e2))
        stop("Matrices must have the same dimensions in order to multiply them.")
    check_valid_matrix(e1)

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

op_str_to_fun <- function(op) {
    return(switch(op, "*"=`*`, "^"=`^`, "/"=`/`, "%%"=`%%`, "%/%"=`%/%`))
}

multiply_csr_by_dvec_elemwise_internal <- function(e1, e2, logical=FALSE,
                                                   X_is_LHS=TRUE, op="*") {
    if (is.matrix(e2)) {
        if (nrow(e1) != nrow(e2) || ncol(e2) != ncol(e2))
            stop("Matrix dimensions do not match. Cannot perform the opertion.")
        if (!logical)
            e2 <- as.numeric(e2)
        else
            e2 <- as.logical(e2)
    }

    if (!NROW(e2)) {
        if (logical)
            return(logical())
        else
            return(numeric())
    }
    if (NROW(e2) > nrow(e1)*ncol(e1)) {
        stop("Vector to multiply with has more entries than matrix.")
    }

    keep_NAs <- !getOption("MatrixExtra.ignore_na", default=FALSE)

    if (!X_is_LHS && keep_NAs && op %in% c("^", "/", "%%", "%/%")) {
        warning("Requested operation is not efficient with a sparse matrix as RHS.")
        e1 <- as(e1, "CsparseMatrix")
        fop <- op_str_to_fun(op)
        return(fop(e2, e1))
    }

    check_valid_matrix(e1)
    take_route_NAs <- (
        !logical && keep_NAs &&
        (   anyNA(e2) ||
            (op %in% c("^", "/", "%%", "%/%") && contains_any_zero(e2)) ||
            (op == "*" && contains_any_inf(e2)) ||
            (op == "^" && contains_any_neg(e2))
        )
    )

    is_coo <- inherits(e1, "TsparseMatrix")
    if (!inherits(e1, c("RsparseMatrix", "TsparseMatrix")))
        throw_internal_error()
    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)


    if (take_route_NAs) {

        if (all(is.na(e2))) {
            warning("Requested operation is not efficient, all values will be NA.")
            e1 <- as(e1, "CsparseMatrix")
            fop <- op_str_to_fun(op)
            if (X_is_LHS)
                return(fop(e1, e2))
            else
                return(fop(e2, e1))
        }

        if (is_coo) {
            e1 <- as.csr.matrix(e1, logical=logical)
            is_coo <- FALSE
        }

        if (inplace_sort)
            e1 <- deepcopy_before_sort(e1)
    }

    if (!logical)
        e2 <- as.numeric(e2)
    else
        e2 <- as.logical(e2)

    if (!inherits(e1, ifelse(logical, "lsparseMatrix", "dsparseMatrix")) ||
        (.hasSlot(e1, "diag") && e1@diag != "N") ||
        inherits(e1, "symmetricMatrix")
    ) {
        if (!is_coo)
            e1 <- as.csr.matrix(e1, logical=logical)
        else
            e1 <- as.coo.matrix(e1, logical=logical)
    }

    if (NROW(e2) == 1L) {

        if (logical) {

            if (is.na(e2)) {
                e1@x <- ifelse(is.na(e1@x), as.logical(NA), FALSE)
                return(e1)
            } else if (!e2) {
                out <- new(ifelse(is_coo, "lgTMatrix", "lgRMatrix"))
                out@Dim <- e1@Dim
                out@Dimnames <- e1@Dimnames
                if (!is_coo)
                    out@p <- rep(0L, nrow(e1)+1L)
                return(out)
            } else {
                return(e1)
            }

        } else if (op == "*") {

            if (is.infinite(e2)) {
                warning("Warning: multiplication by infinite is not efficient for sparse inputs.")
                e1 <- as.csc.matrix(e1)
                return(e1 * e2)
            }

            e1@x <- e1@x * e2

        } else if (op %in% c("/", "%%", "%/%")) {

            fop <- op_str_to_fun(op)
            if (e2 == 0) {
                if (X_is_LHS)
                    warning("Warning: division by zero.")
                if (keep_NAs) {
                    warning("Warning: operation is not efficient for sparse inputs.")
                    e1 <- as.csc.matrix(e1)
                    return(fop(e1, 0))
                } else {
                    if (X_is_LHS)
                        e1@x <- fop(e1@x, e2)
                    else
                        e1@x <- fop(e2, e1@x)
                }
            } else {
                if (X_is_LHS)
                    e1@x <- fop(e1@x, e2)
                else
                    e1@x <- fop(e2, e1@x)
            }

        } else if (op == "^") {

            if (X_is_LHS && keep_NAs && e2 <= 0) {
                warning("Warning: operation is not efficient for sparse inputs.")
                e1 <- as.csc.matrix(e1)
                return(e1 ^ e2)
            }

            if (!X_is_LHS && (e2 == 1 || e2 == 0)) {
                warning("Warning: operation is not efficient for sparse inputs.")
                e1 <- as.csc.matrix(e1)
                return(e2 ^ e1)
            }

            if (X_is_LHS)
                e1@x <- e1@x ^ e2
            else
                e1@x <- e2 ^ e1@x

        } else {

            throw_internal_error()

        }

        return(e1)

    }

    if (take_route_NAs) {
        e1 <- sort_sparse_indices(e1, copy=!inplace_sort)
    }

    if ((nrow(e1) %% length(e2)) != 0)
        warning("Number of elements in vector is not a multiple of matrix dimension.")

    if (take_route_NAs) {
        res <- multiply_csr_by_dvec_with_NAs(
                e1@p, e1@j, e1@x, e2, ncol(e1),
                op=="*", op=="^", op=="/", op=="%%", op=="%/%", X_is_LHS
            )
        X_attr <- attributes(e1)
        X_attr$p <- res$indptr
        X_attr$j <- res$indices
        X_attr$x <- res$values
        attributes(e1) <- X_attr
        return(e1)
    } else {
        if (!logical) {

            if (!is_coo) {
                e1@x <- multiply_csr_by_dvec_no_NAs_numeric(
                    e1@p, e1@j, e1@x, e2, ncol(e1),
                    op=="*", op=="^", op=="/", op=="%%", op=="%/%", X_is_LHS
                )
            } else {
                e1@x <- multiply_coo_by_dense_ignore_NAs_numeric(
                    e1@i, e1@j, e1@x, e2, nrow(e1), ncol(e1),
                    op=="*", op=="^", op=="/", op=="%%", op=="%/%", X_is_LHS
                )
            }

        } else {
            if (!is_coo)
                e1@x <- logicaland_csr_by_dvec_internal(e1@p, e1@j, e1@x, e2, ncol(e1))
            else
                e1@x <- multiply_coo_by_dense_ignore_NAs_logical(e1@i, e1@j, e1@x, e2, nrow(e1), ncol(e1))
        }
        return(e1)
    }
}

multiply_csr_by_dvec_elemwise <- function(e1, e2) {
    if (inherits(e1, "sparseMatrix")) {
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, FALSE))
    } else {
        return(multiply_csr_by_dvec_elemwise_internal(e2, e1, FALSE))
    }
}

logicaland_csr_by_dvec <- function(e1, e2) {
    if (inherits(e1, "sparseMatrix")) {
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, TRUE))
    } else {
        return(multiply_csr_by_dvec_elemwise_internal(e2, e1, TRUE))
    }
}

powerto_csr_by_dvec_elemwise <- function(e1, e2) {
    if (inherits(e1, "sparseMatrix")) {
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, FALSE, op="^", X_is_LHS=TRUE))
    } else {
        return(multiply_csr_by_dvec_elemwise_internal(e2, e1, FALSE, op="^", X_is_LHS=FALSE))
    }
}

divide_csr_by_dvec_elemwise <- function(e1, e2) {
    if (inherits(e1, "sparseMatrix")) {
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, FALSE, op="/", X_is_LHS=TRUE))
    } else {
        return(multiply_csr_by_dvec_elemwise_internal(e2, e1, FALSE, op="/", X_is_LHS=FALSE))
    }
}

intdivide_csr_by_dvec_elemwise <- function(e1, e2) {
    if (inherits(e1, "sparseMatrix")) {
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, FALSE, op="%/%", X_is_LHS=TRUE))
    } else {
        return(multiply_csr_by_dvec_elemwise_internal(e2, e1, FALSE, op="%/%", X_is_LHS=FALSE))
    }
}

divrest_csr_by_dvec_elemwise <- function(e1, e2) {
    if (inherits(e1, "sparseMatrix")) {
        return(multiply_csr_by_dvec_elemwise_internal(e1, e2, FALSE, op="%%", X_is_LHS=TRUE))
    } else {
        return(multiply_csr_by_dvec_elemwise_internal(e2, e1, FALSE, op="%%", X_is_LHS=FALSE))
    }
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="integer"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="numeric"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="logical"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="integer", e2="RsparseMatrix"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="numeric", e2="RsparseMatrix"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="logical", e2="RsparseMatrix"), multiply_csr_by_dvec_elemwise)


#' @rdname operators
#' @export
setMethod("&", signature(e1="RsparseMatrix", e2="integer"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="RsparseMatrix", e2="numeric"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="RsparseMatrix", e2="logical"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="integer", e2="RsparseMatrix"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="numeric", e2="RsparseMatrix"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="logical", e2="RsparseMatrix"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("/", signature(e1="RsparseMatrix", e2="integer"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="RsparseMatrix", e2="numeric"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="RsparseMatrix", e2="logical"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="RsparseMatrix", e2="matrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="integer", e2="RsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="numeric", e2="RsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="logical", e2="RsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="matrix", e2="RsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="RsparseMatrix", e2="integer"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="RsparseMatrix", e2="numeric"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="RsparseMatrix", e2="logical"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="RsparseMatrix", e2="matrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="integer", e2="RsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="numeric", e2="RsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="logical", e2="RsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="matrix", e2="RsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="RsparseMatrix", e2="integer"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="RsparseMatrix", e2="numeric"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="RsparseMatrix", e2="logical"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="RsparseMatrix", e2="matrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="integer", e2="RsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="numeric", e2="RsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="logical", e2="RsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="matrix", e2="RsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="RsparseMatrix", e2="integer"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="RsparseMatrix", e2="numeric"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="RsparseMatrix", e2="logical"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="RsparseMatrix", e2="matrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="integer", e2="RsparseMatrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="numeric", e2="RsparseMatrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="logical", e2="RsparseMatrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="matrix", e2="RsparseMatrix"), powerto_csr_by_dvec_elemwise)


#### Now for COO

#' @rdname operators
#' @export
setMethod("*", signature(e1="TsparseMatrix", e2="integer"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="TsparseMatrix", e2="numeric"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="TsparseMatrix", e2="logical"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="integer", e2="TsparseMatrix"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="numeric", e2="TsparseMatrix"), multiply_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="logical", e2="TsparseMatrix"), multiply_csr_by_dvec_elemwise)


#' @rdname operators
#' @export
setMethod("&", signature(e1="TsparseMatrix", e2="integer"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="TsparseMatrix", e2="numeric"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="TsparseMatrix", e2="logical"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="integer", e2="TsparseMatrix"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="numeric", e2="TsparseMatrix"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("&", signature(e1="logical", e2="TsparseMatrix"), logicaland_csr_by_dvec)

#' @rdname operators
#' @export
setMethod("/", signature(e1="TsparseMatrix", e2="integer"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="TsparseMatrix", e2="numeric"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="TsparseMatrix", e2="logical"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="TsparseMatrix", e2="matrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="integer", e2="TsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="numeric", e2="TsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="logical", e2="TsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("/", signature(e1="matrix", e2="TsparseMatrix"), divide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="TsparseMatrix", e2="integer"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="TsparseMatrix", e2="numeric"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="TsparseMatrix", e2="logical"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="TsparseMatrix", e2="matrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="integer", e2="TsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="numeric", e2="TsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="logical", e2="TsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%%", signature(e1="matrix", e2="TsparseMatrix"), divrest_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="TsparseMatrix", e2="integer"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="TsparseMatrix", e2="numeric"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="TsparseMatrix", e2="logical"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="TsparseMatrix", e2="matrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="integer", e2="TsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="numeric", e2="TsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="logical", e2="TsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("%/%", signature(e1="matrix", e2="TsparseMatrix"), intdivide_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="TsparseMatrix", e2="integer"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="TsparseMatrix", e2="numeric"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="TsparseMatrix", e2="logical"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="TsparseMatrix", e2="matrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="integer", e2="TsparseMatrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="numeric", e2="TsparseMatrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="logical", e2="TsparseMatrix"), powerto_csr_by_dvec_elemwise)

#' @rdname operators
#' @export
setMethod("^", signature(e1="matrix", e2="TsparseMatrix"), powerto_csr_by_dvec_elemwise)

#######



multiply_csr_by_svec_elemwise_internal <- function(X, v) {
    if (!length(v))
        return(numeric())

    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)

    if (length(v) == 1 || length(v) == length(v@i)) {
        if (length(v) == 1) {
            v <- as.numeric(v)
        } else if (.hasSlot(v, "x")) {
            if (inplace_sort)
                v <- deepcopy_before_sort(v)
            v <- sort_sparse_indices(v, copy=!inplace_sort)
            v <- v@x
        } else {
            return(X)
        }
        return(multiply_csr_by_dvec_elemwise(X, v))
    }

    if ((length(v) < nrow(X) && (nrow(X) %% length(v)) != 0) ||
        length(v) > nrow(X)
    ) {
        X <- as.csc.matrix(X, logical=inherits(X, "lsparseMatrix"), binary=inherits(X, "nsparseMatrix"))
        return(X * v)
    }

    check_valid_matrix(X)

    
    if (inplace_sort) {
        X <- deepcopy_before_sort(X)
        v <- deepcopy_before_sort(v)
    }
    X <- as.csr.matrix(X)
    X <- sort_sparse_indices(X, copy=!inplace_sort)

    v <- as.sparse.vector(v, binary=inherits(v, "nsparseVector"))
    v <- sort_sparse_indices(v, copy=!inplace_sort)

    keep_NAs <- !getOption("MatrixExtra.ignore_na", default=FALSE)

    if (.hasSlot(v, "x")) {
        if (keep_NAs) {
            res <- multiply_csr_by_svec_keep_NAs(X@p, X@j, X@x, v@i, v@x, ncol(X), length(v))
        } else {
            res <- multiply_csr_by_svec_no_NAs(X@p, X@j, X@x, v@i, v@x, length(v))
        }
    } else {
        res <- multiply_csr_by_svec_no_NAs(X@p, X@j, X@x, v@i, numeric(), length(v))
    }
    out <- new("dgRMatrix")
    out@Dim <- X@Dim
    out@Dimnames <- X@Dimnames
    out@p <- res$indptr
    out@j <- res$indices
    out@x <- res$values
    return(out)
}

multiply_csr_by_svec_elemwise <- function(e1, e2) {
    if (inherits(e1, "sparseMatrix")) {
        return(multiply_csr_by_svec_elemwise_internal(e1, e2))
    } else {
        return(multiply_csr_by_svec_elemwise_internal(e2, e1))
    }
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="RsparseMatrix", e2="sparseVector"), multiply_csr_by_svec_elemwise)

#' @rdname operators
#' @export
setMethod("*", signature(e1="sparseVector", e2="RsparseMatrix"), multiply_csr_by_svec_elemwise)


multiply_elemwise_dense_by_svec_internal <- function(e1, e2) {
    if (!NROW(e1) || !NCOL(e2) || !length(e2))
        return(matrix())

    keep_NAs <- !getOption("MatrixExtra.ignore_na", default=FALSE)
    inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)
    if (keep_NAs && inplace_sort)
        e2 <- deepcopy_before_sort(e2)
    e2 <- as.sparse.vector(e2)
    if (keep_NAs)
        e2 <- sort_sparse_indices(e2, copy=!inplace_sort)

    if (typeof(e1) == "double") {
        res <- multiply_elemwise_dense_by_svec_numeric(e1, e2@i, e2@x, length(e2), keep_NAs)
    } else if (typeof(e1) == "integer") {
        res <- multiply_elemwise_dense_by_svec_integer(e1, e2@i, e2@x, length(e2), keep_NAs)
    } else if (typeof(e1) == "logical") {
        res <- multiply_elemwise_dense_by_svec_logical(e1, e2@i, e2@x, length(e2), keep_NAs)
    } else if (inherits(e1, "float32")) {
        ### TODO: revisit this
        if (is.vector(e1@Data))
            e1@Data <- as.matrix(e1@Data)
        res <- multiply_elemwise_dense_by_svec_float32(e1@Data, e2@i, e2@x, length(e2), keep_NAs)
    } else {
        mode(e1) <- "double"
        return(multiply_elemwise_dense_by_svec_internal(e1, e2))
    }

    if ("X_dense" %in% names(res))
        return(res$X_dense)

    out <- new("dgRMatrix")
    out@Dim <- dim(e1)
    if (!is.null(rownames(e1)))
        rownames(out) <- rownames(e1)
    if (!is.null(colnames(e1)))
        colnames(out) <- colnames(e1)
    out@p <- res$indptr
    out@j <- res$indices
    out@x <- res$values
    return(out)
}

multiply_elemwise_dense_by_svec <- function(e1, e2) {
    if (is.matrix(e1))
        return(multiply_elemwise_dense_by_svec_internal(e1, e2))
    else
        return(multiply_elemwise_dense_by_svec_internal(e2, e1))
}

#' @rdname operators
#' @export
setMethod("*", signature(e1="matrix", e2="sparseVector"), multiply_elemwise_dense_by_svec)

#' @rdname operators
#' @export
setMethod("*", signature(e1="sparseVector", e2="matrix"), multiply_elemwise_dense_by_svec)

#' @rdname operators
#' @export
setMethod("*", signature(e1="float32", e2="sparseVector"), multiply_elemwise_dense_by_svec)

#' @rdname operators
#' @export
setMethod("*", signature(e1="sparseVector", e2="float32"), multiply_elemwise_dense_by_svec)

