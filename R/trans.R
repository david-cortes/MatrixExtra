t_csc_to_csr <- function(X) {
    X_attr <- attributes(X)
    X_attr$class <- gsub("CMatrix", "RMatrix", X_attr$class, ignore.case=FALSE)
    X_attr$j <- X_attr$i
    X_attr$i <- NULL
    X_attr$Dim <- rev(X_attr$Dim)
    X_attr$Dimnames <- rev(X_attr$Dimnames)
    if ("factors" %in% names(X_attr))
        X_attr$factors <- list()
    if ("uplo" %in% names(X_attr))
        X_attr$uplo <- ifelse(X_attr$uplo == "L", "U", "L")
    attributes(X) <- X_attr
    return(X)
}

t_csr_to_csc <- function(X) {
    X_attr <- attributes(X)
    X_attr$class <- gsub("RMatrix", "CMatrix", X_attr$class, ignore.case=FALSE)
    X_attr$i <- X_attr$j
    X_attr$j <- NULL
    X_attr$Dim <- rev(X_attr$Dim)
    X_attr$Dimnames <- rev(X_attr$Dimnames)
    if ("factors" %in% names(X_attr))
        X_attr$factors <- list()
    if ("uplo" %in% names(X_attr))
        X_attr$uplo <- ifelse(X_attr$uplo == "L", "U", "L")
    attributes(X) <- X_attr
    return(X)
}

t_coo_to_coo <- function(X) {
    X_attr <- attributes(X)
    temp <- X_attr$i
    X_attr$i <- X_attr$j
    X_attr$j <- temp
    X_attr$Dim <- rev(X_attr$Dim)
    X_attr$Dimnames <- rev(X_attr$Dimnames)
    if ("factors" %in% names(X_attr))
        X_attr$factors <- list()
    if ("uplo" %in% names(X_attr))
        X_attr$uplo <- ifelse(X_attr$uplo == "L", "U", "L")
    attributes(X) <- X_attr
    return(X)
}

t_deep_internal <- function(x) {
    check_valid_matrix(x)
    orig_class <- class(x)
    x <- as(x, "TsparseMatrix")
    x <- t_shallow(x)
    if (grepl("CMatrix", orig_class))
        x <- as(x, "CsparseMatrix")
    else
        x <- as(x, "RsparseMatrix")
    return(x)
}

t_masked_csc <- function(x) {
    if (getOption("MatrixExtra.fast_transpose", default=FALSE)) {
        return(t_shallow(x))
    } else {
        return(t_deep_internal(x))
    }
}

t_masked_csr <- function(x) {
    if (getOption("MatrixExtra.fast_transpose", default=FALSE)) {
        return(t_shallow(x))
    } else {
        return(t_deep_internal(x))
    }
}

t_masked_coo <- function(x) {
    return(t_shallow(x))
}

#' @title Transpose a sparse matrix by changing its format
#' @description Transposes a sparse matrix in CSC (a.k.a. "CsparseMatrix")
#' or CSR (a.k.a. "RsparseMatrix") formats by converting it to the opposite format
#' (i.e. CSC -> CSR, CSR -> CSC).
#'
#' This implies only a shallow copy (i.e. it's much faster), as the only necessary thing to make
#' such transpose operation is to swap the number of rows and columns and change the class
#' of the object (all data remains the same), avoiding any deep copying and
#' format conversion as when e.g. creating a CSC transpose of a CSC matrix.
#'
#' If the input is neither a CSR not CSC matrix, it will just call the generic `t()` method.
#' 
#' Also provided is a function `t_deep` which outputs a transpose with the same storage order.
#' @details \bold{Important:} When loading this package (`library(MatrixExtra)`), it will
#' change the behavior of `t(sparseMatrix)` towards calling `t_shallow`.
#' 
#' This makes it more efficient, but has the potential of breaking existing code in other
#' packages, particularly in the `Matrix` package itself when calling some arbitrary
#' function or method which would internally transpose a CSC matrix and rely on the assumption
#' that its output is also CSC.
#' 
#' This behavior can be changed through \link{restore_old_matrix_behavior} or
#' the package options (e.g. `options("MatrixExtra.fast_transpose" = FALSE)` -
#' ee \link{MatrixExtra-options}) to have `t_deep` as the default, just like in `Matrix`.
#' 
#' Additionally, under the new behavior (`t_shallow` as the default for `t`),
#' transposing a `sparseVector` object will yield a CSR matrix ("RsparseMatrix"),
#' which differs from `Matrix` that would yield a COO matrix ("TsparseMatrix").
#' @param x A sparse matrix. If `x` is of a different type, will just invoke its generic
#' `t()` method.
#' @returns The transpose of `x` (rows become columns and columns become rows),
#' but in the opposite format (CSC -> CSR, CSR -> CSC); or the same format if calling `t_deep`.
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(3, 4, .5, repr="C")
#' inherits(X, "CsparseMatrix")
#' Xtrans <- t_shallow(X)
#' inherits(Xtrans, "RsparseMatrix")
#' nrow(X) == ncol(Xtrans)
#' ncol(X) == nrow(Xtrans)
#'
#' Xorig <- t_shallow(Xtrans)
#' inherits(Xorig, "CsparseMatrix")
#' inherits(t_deep(Xtrans), "RsparseMatrix")
#' 
#' ### Important!!!
#' ### This package makes 't_shallow' the default
#' set_new_matrix_behavior()
#' inherits(X, "CsparseMatrix")
#' inherits(t(X), "RsparseMatrix")
#' 
#' ### Can be changed back to 't_deep' like this:
#' restore_old_matrix_behavior()
#' inherits(t(X), "CsparseMatrix")
#' @export
t_shallow <- function(x) {

    if (inherits(x, "sparseMatrix")) {
        if (inherits(x, "CsparseMatrix")) {
            return(t_csc_to_csr(x))
        } else if (inherits(x, "RsparseMatrix")) {
            return(t_csr_to_csc(x))
        } else if (inherits(x, "TsparseMatrix")) {
            return(t_coo_to_coo(x))
        } else {
            throw_internal_error()
        }
    }
    return(t(x))
}

#' @rdname t_shallow
#' @export
t_deep <- function(x) {
    if (inherits(x, c("RsparseMatrix", "CsparseMatrix")))
        return(t_deep_internal(x))
    else
        return(t(x))
}

#' @rdname t_shallow
#' @export
setMethod("t", signature(x="RsparseMatrix"), t_masked_csr)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="CsparseMatrix"), t_masked_csc)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="TsparseMatrix"), t_masked_coo)

#' @rdname t_shallow
#' @export
setMethod("t", signature(x="dgCMatrix"), t_masked_csc)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="ngCMatrix"), t_masked_csc)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="lgCMatrix"), t_masked_csc)

#' @rdname t_shallow
#' @export
setMethod("t", signature(x="dtCMatrix"), t_masked_csc)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="ntCMatrix"), t_masked_csc)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="ltCMatrix"), t_masked_csc)

#' @rdname t_shallow
#' @export
setMethod("t", signature(x="dsCMatrix"), t_masked_csc)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="nsCMatrix"), t_masked_csc)
#' @rdname t_shallow
#' @export
setMethod("t", signature(x="lsCMatrix"), t_masked_csc)

#' @rdname t_shallow
#' @export
setMethod("t", signature(x="sparseVector"), function(x) {
    if (getOption("MatrixExtra.fast_transpose", default=FALSE))
        return(as.csr.matrix(x, binary=inherits(x, "nsparseVector"), logical=inherits(x, "lsparseVector")))
    else
        return(as.coo.matrix(x, binary=inherits(x, "nsparseVector"), logical=inherits(x, "lsparseVector")))
})
