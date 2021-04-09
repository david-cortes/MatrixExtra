#' @name csr-linalg
#' @title Linear Algebra functions for CSR matrices
#' @description Short wrappers around some linear algebra operators from `Matrix`
#' that take CSC matrices, adapted to work for CSR matrices without involving any
#' data duplication or deep format conversion, thus saving time and memory.
#' @param x A sparse matrix in CSR format.
#' @param type Type of the norm to calculate (see \link[Matrix]{norm}).
#' @param value Replacement value for the matrix diagonal.
#' @param ... Extra arguments to pass to `norm`
#' @return The same value that `Matrix` would return for CSC matrices.
NULL

norm_csr <- function(x, type="O", ...) {
    if (missing(type))
        type <- "O"
    if (typeof(type) != "character")
        stop("'type' must be a single character variable (see ?norm).")
    type <- type[1L]
    allowed_types <- c("O", "o", "1",  "I", "i", "F", "f", "M", "m", "2")
    if (!(type %in% allowed_types))
        stop(sprintf("Invalid norm type. Allowed values: %s", paste(allowed_types, sep=", ")))

    if (type %in% c("O", "o", "1")) {
        return(norm(t_shallow(x), "I", ...))
    } else if (type %in% c("I", "i")) {
        ### TODO: this one can be done more efficiently
        return(norm(t_shallow(x), "O", ...))
    } else {
        return(norm(t_shallow(x), type, ...))
    }
}

#' @rdname csr-linalg
#' @export
setMethod("norm", signature(x="RsparseMatrix", type="character"), norm_csr)

#' @rdname csr-linalg
#' @export
setMethod("norm", signature(x="RsparseMatrix", type="missing"), norm_csr)

diag_csr <- function(x) {
    return(diag(t_shallow(x)))
}

assign_diag_csr <- function(x, value) {
    x <- t_shallow(x)
    diag(x) <- value
    return(t_shallow(x))
}

#' @rdname csr-linalg
#' @export
setMethod("diag", signature(x="RsparseMatrix"), diag_csr)

#' @rdname csr-linalg
#' @export
setMethod("diag<-", signature(x="RsparseMatrix", value="ANY"), assign_diag_csr)
