#' @name mathematical-functions
#' @title Mathematical functions for CSR and COO matrices
#' @description Implements some mathematical functions for CSR
#' (a.k.a. "RsparseMatrix") and COO (a.k.a. "TsparseMatrix") matrix types
#' without converting them to CSC matrices in the process.
#'
#' The functions here reduce to applying the same operation over the non-zero elements only,
#' and as such do not benefit from any storage format conversion as done implicitly
#' in the `Matrix` package.
#' @param x A CSR or COO matrix.
#' @param digits See \link{round} and \link{signif}. If passing more than one value,
#' will call the corresponding function from the `Matrix` package, which implies first
#' converting `x` to CSC format.
#' @return A CSR or COO matrix depending on the input. They will be of the `dg` type
#' (`dgRMatrix` or `dgTMatrix`).
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- as.csr.matrix(rsparsematrix(4, 3, .4))
#' abs(X)
#' sqrt(X^2)
#' ### This will output CSC
#' round(X, 1:2)
NULL

#' @rdname mathematical-functions
#' @export
setMethod("sqrt", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "U")) {
        x <- as.csr.matrix(x)
    }
    x@x <- sqrt(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("sqrt", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "U")) {
        x <- as.coo.matrix(x)
    }
    x@x <- sqrt(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("abs", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix"))
        x <- as.csr.matrix(x)
    x@x <- abs(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("abs", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix"))
        x <- as.coo.matrix(x)
    x@x <- abs(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("log1p", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- log1p(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("log1p", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- log1p(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("sin", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- sin(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("sin", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- sin(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("tan", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- tan(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("tan", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- tan(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("tanh", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- tanh(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("tanh", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- tanh(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("tanpi", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix"))
        x <- as.csr.matrix(x)
    x@x <- tanpi(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("tanpi", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix"))
        x <- as.coo.matrix(x)
    x@x <- tanpi(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("sinh", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- sinh(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("sinh", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- sinh(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("atanh", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- atanh(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("atanh", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- atanh(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("expm1", signature(x="RsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- expm1(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("expm1", signature(x="TsparseMatrix"), function(x) {
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- expm1(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("sign", signature(x="RsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.csr.matrix(x))
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x)
    x@x <- sign(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("sign", signature(x="TsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.coo.matrix(x))
    if (!inherits(x, "dsparseMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x)
    x@x <- sign(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("ceiling", signature(x="RsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.csr.matrix(x))
    if (!inherits(x, "dsparseMatrix"))
        x <- as.csr.matrix(x)
    x@x <- ceiling(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("ceiling", signature(x="TsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.coo.matrix(x))
    if (!inherits(x, "dsparseMatrix"))
        x <- as.coo.matrix(x)
    x@x <- ceiling(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("floor", signature(x="RsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.csr.matrix(x))
    if (!inherits(x, "dsparseMatrix"))
        x <- as.csr.matrix(x)
    x@x <- floor(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("floor", signature(x="TsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.coo.matrix(x))
    if (!inherits(x, "dsparseMatrix"))
        x <- as.coo.matrix(x)
    x@x <- floor(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("trunc", signature(x="RsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.csr.matrix(x))
    if (!inherits(x, "dsparseMatrix"))
        x <- as.csr.matrix(x)
    x@x <- trunc(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("trunc", signature(x="TsparseMatrix"), function(x) {
    if (inherits(x, "nsparseMatrix"))
        return(as.coo.matrix(x))
    if (!inherits(x, "dsparseMatrix"))
        x <- as.coo.matrix(x)
    x@x <- trunc(x@x)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("round", signature(x="RsparseMatrix", digits="ANY"), function(x, digits) {
    if (inherits(x, "nsparseMatrix"))
        return(as.csr.matrix(x))
    if (!missing(digits) && NROW(digits) != 1L) {
        x <- as(x, "CsparseMatrix")
        return(round(x, digits))
    }
    if (!inherits(x, "dsparseMatrix"))
        x <- as.csr.matrix(x)
    x@x <- round(x@x, digits)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("round", signature(x="TsparseMatrix", digits="ANY"), function(x, digits) {
    if (inherits(x, "nsparseMatrix"))
        return(as.coo.matrix(x))
    if (!missing(digits) && NROW(digits) != 1L) {
        x <- as(x, "CsparseMatrix")
        return(round(x, digits))
    }
    if (!inherits(x, "dsparseMatrix"))
        x <- as.coo.matrix(x)
    x@x <- round(x@x, digits)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("signif", signature(x="RsparseMatrix", digits="ANY"), function(x, digits) {
    if (inherits(x, "nsparseMatrix"))
        return(as.csr.matrix(x))
    if (!missing(digits) && NROW(digits) != 1L) {
        x <- as(x, "CsparseMatrix")
        return(signif(x, digits))
    }
    if (!inherits(x, "dsparseMatrix"))
        x <- as.csr.matrix(x)
    x@x <- signif(x@x, digits)
    return(x)
})

#' @rdname mathematical-functions
#' @export
setMethod("signif", signature(x="TsparseMatrix", digits="ANY"), function(x, digits) {
    if (inherits(x, "nsparseMatrix"))
        return(as.coo.matrix(x))
    if (!missing(digits) && NROW(digits) != 1L) {
        x <- as(x, "CsparseMatrix")
        return(signif(x, digits))
    }
    if (!inherits(x, "dsparseMatrix"))
        x <- as.coo.matrix(x)
    x@x <- signif(x@x, digits)
    return(x)
})
