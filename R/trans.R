t_csc_to_csr <- function(X) {
    out <- new(gsub("CMatrix", "RMatrix", class(X)[1L], ignore.case=FALSE))
    out@Dim <- rev(X@Dim)
    out@Dimnames <- rev(X@Dimnames)
    out@p <- X@p
    out@j <- X@i
    if (.hasSlot(X, "x"))
        out@x <- X@x
    if (.hasSlot(X, "diag"))
        out@diag <- X@diag
    if (.hasSlot(X, "uplo"))
        out@uplo <- ifelse(X@uplo == "L", "U", "L")
    return(out)
}

t_csr_to_csc <- function(X) {
    out <- new(gsub("RMatrix", "CMatrix", class(X)[1L], ignore.case=FALSE))
    out@Dim <- rev(X@Dim)
    out@Dimnames <- rev(X@Dimnames)
    out@p <- X@p
    out@i <- X@j
    if (.hasSlot(X, "x"))
        out@x <- X@x
    if (.hasSlot(X, "diag"))
        out@diag <- X@diag
    if (.hasSlot(X, "uplo"))
        out@uplo <- ifelse(X@uplo == "L", "U", "L")
    return(out)
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
