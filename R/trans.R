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

#' @title Transpose a sparse matrix by changing its format
#' @description Transposes a sparse matrix in CSC (a.k.a. "CsparseMatrix")
#' or CSR (a.k.a. "RsparseMatrix") formats by converting it to the opposite format
#' (i.e. CSC -> CSR, CSR -> CSC).
#'
#' This implies only a shallow copy (i.e. it's faster), as the only necessary thing to make
#' such transpose operation is to swap the number of rows and columns and change the class
#' of the object (all other slots remain the same), avoiding any deep copying and
#' format conversion as when e.g. creating a CSC transpose of a CSC matrix.
#'
#' If the input is neither a CSR not CSC matrix, it will just call the generic `t()` method.
#' @param X A sparse matrix in CSC (`dgCMatrix` or `ngCMatrix`) or CSR (`dgRMatrix` or `ngRMatrix`) formats. If `X` is of a different type, will just invoke its generic
#' `t()` method.
#' @returns The transpose of `X` (rows become columns and columns become rows),
#' but in the opposite format (CSC -> CSR, CSR -> CSC).
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(3, 4, .5)
#' inherits(X, "CsparseMatrix")
#' Xtrans <- t_shallow(X)
#' inherits(Xtrans, "RsparseMatrix")
#' nrow(X) == ncol(Xtrans)
#' ncol(X) == nrow(Xtrans)
#'
#' Xorig <- t_shallow(Xtrans)
#' inherits(Xorig, "CsparseMatrix")
#' @export
t_shallow <- function(X) {
    if (inherits(X, "CsparseMatrix")) {
        return(t_csc_to_csr(X))
    } else if (inherits(X, "RsparseMatrix")) {
        return(t_csr_to_csc(X))
    } else {
        return(t(X))
    }
}
