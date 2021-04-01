#' @name cbind2-method
#' @title Concatenate COO matrices by columns
#' @description `cbind2` method for the sparse COO matrix and sparse vector classes from `Matrix`,
#' taking the most efficient route for the concatenation according to the input types.
#' @param x First matrix to concatenate.
#' @param y Second matrix to concatenate.
#' @return A COO matrix.
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(3, 4, .3)
#' X <- as(X, "TsparseMatrix")
#' inherits(cbind2(X, X), "dgTMatrix")
NULL

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="TsparseMatrix", y="TsparseMatrix"), function(x, y) {
    return(t(rbind2_tri(t(x), t(y))))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="TsparseMatrix", y="sparseVector"), function(x, y) {
    return(t(rbind2_tri_vec(t(x), y, TRUE)))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="sparseVector", y="TsparseMatrix"), function(x, y) {
    return(t(rbind2_tri_vec(t(y), x, FALSE)))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="CsparseMatrix", y="sparseVector"), function(x, y) {
    return(t_shallow(rbind2_generic(t_shallow(x), y)))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="sparseVector", y="CsparseMatrix"), function(x, y) {
    return(t_shallow(rbind2_generic(x, t_shallow(y))))
})

#' @rdname cbind2-method
#' @export
setMethod("cbind2", signature(x="sparseVector", y="sparseVector"), function(x, y) {
    return(t_shallow(rbind2_generic(x, y)))
})
