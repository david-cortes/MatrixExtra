#' @name conversions
#' @title Conversions between matrix types
#' @description Convenience functions for converting to different sparse matrix formats,
#' between pairs of classes which might not be supported in the `Matrix` package.
#'
#' These come in the form of explicit functions 'as.<type>.matrix' (see below),
#' as well as registered conversion methods to use along with `as(object, type)`, adding
#' extra conversion routes which are missing in the `Matrix` package for output
#' types `dgRMatrix`, `lgRMatrix`, and `ngRMatrix`.
#' @details The functions internally might use some routes of `as(x, "?sparseMatrix")`, so they might work
#' with other object classes if they register a conversion method for `Matrix` base
#' types.
#'
#' When passed a vector, the functions `as.csr.matrix` and `as.coo.matrix` will
#' assume that it is a row vector, while `as.csc.matrix` will assume it's a column vector.
#' @param x A matrix which is to be converted to a different format.
#'
#' Supported input types are:\itemize{
#' \item Sparse matrices from `Matrix` package, in any format.
#' \item Sparse vectors from `Matrix` in any format.
#' \item Dense matrices from base R (class `matrix`).
#' \item Dense vectors from base R (classes `numeric`, `integer`, `logical`).
#' \item Dense matrix or vector from package `float` (class `float32`).
#' \item `data.frame`, `data.table`, and `tibble`.
#' }
#' @param binary Whether the result should be a binary-only matrix/vector (inheriting from
#' class `nsparseMatrix`/`nsparseVector` - these don't have slot `x`).
#' Can only pass one of `binary` or `logical`.
#' @param logical Whether the result should be a matrix/vector with logical (boolean) type
#' (inheriting from `lsparseMatrix`/`lsparseVector`).
#' Can only pass one of `binary` or `logical`.
#' @param integer Whether the result should be a vector with integer type ('isparseVector').
#' @param sort Whether to sort the indices in case they are not sorted. Note that it will
#' perform deep copies of the indices and values along the way.
#' @return A sparse matrix/vector, with format:\itemize{
#' \item CSR (a.k.a. `RsparseMatrix`) when calling `as.csr.matrix`
#' (class `dgRMatrix`, `ngRMatrix`, or `lgRMatrix`, depending on parameters `binary` and `logical`).
#' \item CSC (a.k.a. `CsparseMatrix`) when calling `as.csc.matrix`
#' (class `dgCMatrix`, `ngCMatrix`, or `lgCMatrix`, depending on parameters `binary` and `logical`).
#' \item COO (a.k.a. `TsparseMatrix`) when calling `as.coo.matrix`
#' (class `dgTMatrix`, `ngTMatrix`, or `lgTMatrix`, depending on parameters `binary` and `logical`).
#' \item sparse vector (class dependant on input) when calling `as.sparse.vector`.
#' }
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#'
#' m.coo <- as(matrix(1:3), "TsparseMatrix")
#' as.csr.matrix(m.coo)
#' as.csr.matrix(1:3) # <- assumes it's a row vector
#' as.csc.matrix(1:3) # <- assumes it's a column vector
#'
#' ### Using the new conversion methods
#' ### (these would fail if 'MatrixExtra' is not loaded)
#' as(matrix(1:3), "ngRMatrix")
#' as(as.csc.matrix(m.coo), "dgRMatrix")
NULL

#' @rdname conversions
#' @export
as.csr.matrix <- function(x, binary=FALSE, logical=FALSE, sort=FALSE) {
    if (binary && logical)
        stop("Can pass only one of 'binary' or 'logical'.")

    if ((inherits(x, "dgRMatrix") && !binary && !logical) ||
        (inherits(x, "ngRMatrix") && binary) ||
        (inherits(x, "lgRMatrix") && logical)) {
        return(x)
    }

    if (inherits(x, "float32"))
        x <- float::dbl(x)

    if (inherits(x, "data.frame"))
        x <- as.matrix(x)

    if (inherits(x, c("numeric", "integer", "logical")))
        x <- matrix(x, nrow=1L)


    if (!binary && !logical) {
        target_class <- "dgRMatrix"
    } else if (binary) {
        target_class <- "ngRMatrix"
    } else {
        target_class <- "lgRMatrix"
    }

    if (inherits(x, "sparseVector")) {
        X.csr <- new(target_class)
        X.csr@Dim <- as.integer(c(1L, x@length))
        X.csr@p <- c(0L, length(x@i))
        X.csr@j <- as.integer(x@i) - 1L
        if (!binary) {

            if (inherits(x, "dsparseVector")) {
                if (!logical)
                    X.csr@x <- x@x
                else
                    X.csr@x <- as.logical(x@x)
            } else if (inherits(x, "isparseVector")) {
                if (!logical)
                    X.csr@x <- as.numeric(x@x)
                else
                    X.csr@x <- as.logical(x@x)
            } else if (inherits(x, "lsparseVector")) {
                    if (!logical)
                        X.csr@x <- as.numeric(x@x)
                    else
                        X.csr@x <- x@x
            } else {
                if (!logical)
                    X.csr@x <- rep(1., length(x@i))
                else
                    X.csr@x <- rep(TRUE, length(x@i))
            }

        }
        x <- X.csr
    }

    if (!inherits(x, "RsparseMatrix"))
        x <- as(x, "RsparseMatrix")

    if (inherits(x, c("symmetricMatrix", "triangularMatrix"))) {
        x_trans <- t_shallow(x)
        if (!inherits(x_trans, "dsparseMatrix"))
            x_trans <- as(x_trans, "dsparseMatrix")
        x_trans <- as(x_trans, "dgCMatrix")
        x <- t_shallow(x_trans)
    }


    if (!binary && !logical && !inherits(x, "dgRMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "dgRMatrix"
        if (.hasSlot(x, "x"))
            X_attr$x <- as.numeric(X_attr$x)
        else
            X_attr$x <- rep(1., length(X_attr$j))
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (logical && !inherits(x, "lgRMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "lgRMatrix"
        if (.hasSlot(x, "x"))
            X_attr$x <- as.logical(X_attr$x)
        else
            X_attr$x <- rep(TRUE, length(X_attr$j))
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (binary && !inherits(x, "ngRMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "ngRMatrix"
        if ("x" %in% names(X_attr))
            X_attr$x <- NULL
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (sort) X <- sort_sparse_indices(X, copy=TRUE)
    return(x)
}

#' @rdname conversions
#' @export
as.csc.matrix <- function(x, binary=FALSE, logical=FALSE, sort=FALSE) {
    if (binary && logical)
        stop("Can pass only one of 'binary' or 'logical'.")

    if ((inherits(x, "dgCMatrix") && !binary && !logical) ||
        (inherits(x, "ngCMatrix") && binary) ||
        (inherits(x, "lgCMatrix") && logical)) {
        return(x)
    }

    if (inherits(x, "float32"))
        x <- float::dbl(x)

    if (inherits(x, c("numeric", "integer", "logical", "data.frame")))
        x <- as.matrix(x)

    if (!inherits(x, "CsparseMatrix"))
        x <- as(x, "CsparseMatrix")

    if (inherits(x, c("symmetricMatrix", "triangularMatrix"))) {
        if (!inherits(x, "dsparseMatrix"))
            x <- as(x, "dsparseMatrix")
        x <- as(x, "dgCMatrix")
    }

    if (!binary && !logical && !inherits(x, "dgCMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "dgCMatrix"
        if (.hasSlot(x, "x"))
            X_attr$x <- as.numeric(X_attr$x)
        else
            X_attr$x <- rep(1., length(X_attr$i))
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (logical && !inherits(x, "lgCMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "lgCMatrix"
        if (.hasSlot(x, "x"))
            X_attr$x <- as.logical(X_attr$x)
        else
            X_attr$x <- rep(TRUE, length(X_attr$i))
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (binary && !inherits(x, "ngCMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "ngCMatrix"
        if ("x" %in% names(X_attr))
            X_attr$x <- NULL
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (sort) X <- sort_sparse_indices(X, copy=TRUE)
    return(x)
}

#' @rdname conversions
#' @export
as.coo.matrix <- function(x, binary=FALSE, logical=FALSE, sort=FALSE) {
    if (binary && logical)
        stop("Can pass only one of 'binary' or 'logical'.")

    if ((inherits(x, "dgTMatrix") && !binary && !logical) ||
        (inherits(x, "ngTMatrix") && binary) ||
        (inherits(x, "lgTMatrix") && logical)) {
        return(x)
    }

    if (inherits(x, "float32"))
        x <- float::dbl(x)

    if (inherits(x, c("data.frame")))
        x <- as.matrix(x)

    if (inherits(x, c("numeric", "integer", "logical")))
        x <- matrix(x, nrow=1L)

    if (inherits(x, "sparseVector"))
        x <- as.csr.matrix(x)

    if (!inherits(x, "TsparseMatrix"))
        x <- as(x, "TsparseMatrix")

    if (inherits(x, c("symmetricMatrix", "triangularMatrix"))) {
        if (!inherits(x, "dsparseMatrix"))
            x <- as(x, "dsparseMatrix")
        x <- as(x, "dgTMatrix")
    }

    if (!binary && !logical && !inherits(x, "dgTMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "dgTMatrix"
        if (.hasSlot(x, "x"))
            X_attr$x <- as.numeric(X_attr$x)
        else
            X_attr$x <- rep(1., length(X_attr$j))
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (logical && !inherits(x, "lgTMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "lgTMatrix"
        if (.hasSlot(x, "x"))
            X_attr$x <- as.logical(X_attr$x)
        else
            X_attr$x <- rep(TRUE, length(X_attr$j))
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (binary && !inherits(x, "ngTMatrix")) {
        X_attr <- attributes(x)
        X_attr$class <- "ngTMatrix"
        if ("x" %in% names(X_attr))
            X_attr$x <- NULL
        if ("diag" %in% names(X_attr))
            X_attr$diag <- NULL
        if ("uplo" %in% names(X_attr))
            X_attr$uplo <- NULL
        attributes(x) <- X_attr
    }

    if (sort) X <- sort_sparse_indices(X, copy=TRUE)
    return(x)
}

#' @rdname conversions
#' @export
as.sparse.vector <- function(x, binary=FALSE, logical=FALSE, integer=FALSE) {
    if ((binary && logical) || (logical && integer) || (binary && integer))
        stop("Can pass at most one of 'binary', 'logical', 'integer'.")

    if (inherits(x, "data.frame"))
        x <- as.matrix(x)
    if (inherits(x, "float32"))
        x <- float::dbl(x)

    x <- as(x, "sparseVector")

    if (!binary && !logical && !integer) {
        if (!inherits(x, "dsparseVector"))
            x <- as(x, "dsparseVector")
    } else if (integer) {
        if (!inherits(x, "isparseVector"))
            x <- as(x, "isparseVector")
    } else if (logical) {
        if (!inherits(x, "lsparseVector"))
            x <- as(x, "lsparseVector")
    } else if (binary) {
        if (!inherits(x, "nsparseVector"))
            x <- as(x, "nsparseVector")
    }

    return(x)
}

#' @export
setAs("sparseMatrix", "dgRMatrix", function(from) as.csr.matrix(from))
#' @export
setAs("sparseMatrix", "lgRMatrix", function(from) as.csr.matrix(from, logical=TRUE))
#' @export
setAs("sparseMatrix", "ngRMatrix", function(from) as.csr.matrix(from, binary=TRUE))

#' @export
setAs("matrix", "dgRMatrix", function(from) as.csr.matrix(from))
#' @export
setAs("matrix", "lgRMatrix", function(from) as.csr.matrix(from, logical=TRUE))
#' @export
setAs("matrix", "ngRMatrix", function(from) as.csr.matrix(from, binary=TRUE))

#' @export
setAs("sparseVector", "dgRMatrix", function(from) as.csr.matrix(from))
#' @export
setAs("sparseVector", "lgRMatrix", function(from) as.csr.matrix(from, logical=TRUE))
#' @export
setAs("sparseVector", "ngRMatrix", function(from) as.csr.matrix(from, binary=TRUE))

#' @export
setAs("numeric", "dgRMatrix", function(from) as.csr.matrix(from))
#' @export
setAs("integer", "dgRMatrix", function(from) as.csr.matrix(from))
#' @export
setAs("logical", "dgRMatrix", function(from) as.csr.matrix(from))

#' @export
setAs("numeric", "lgRMatrix", function(from) as.csr.matrix(from, logical=TRUE))
#' @export
setAs("integer", "lgRMatrix", function(from) as.csr.matrix(from, logical=TRUE))
#' @export
setAs("logical", "lgRMatrix", function(from) as.csr.matrix(from, logical=TRUE))

#' @export
setAs("numeric", "ngRMatrix", function(from) as.csr.matrix(from, binary=TRUE))
#' @export
setAs("integer", "ngRMatrix", function(from) as.csr.matrix(from, binary=TRUE))
#' @export
setAs("logical", "ngRMatrix", function(from) as.csr.matrix(from, binary=TRUE))



#' @export
setAs("sparseMatrix", "dgCMatrix", function(from) as.csc.matrix(from))
#' @export
setAs("sparseMatrix", "lgCMatrix", function(from) as.csc.matrix(from, logical=TRUE))
#' @export
setAs("sparseMatrix", "ngCMatrix", function(from) as.csc.matrix(from, binary=TRUE))

#' @export
setAs("matrix", "dgCMatrix", function(from) as.csc.matrix(from))
#' @export
setAs("matrix", "lgCMatrix", function(from) as.csc.matrix(from, logical=TRUE))
#' @export
setAs("matrix", "ngCMatrix", function(from) as.csc.matrix(from, binary=TRUE))

#' @export
setAs("sparseVector", "dgCMatrix", function(from) as.csc.matrix(from))
#' @export
setAs("sparseVector", "lgCMatrix", function(from) as.csc.matrix(from, logical=TRUE))
#' @export
setAs("sparseVector", "ngCMatrix", function(from) as.csc.matrix(from, binary=TRUE))

#' @export
setAs("numeric", "dgCMatrix", function(from) as.csc.matrix(from))
#' @export
setAs("integer", "dgCMatrix", function(from) as.csc.matrix(from))
#' @export
setAs("logical", "dgCMatrix", function(from) as.csc.matrix(from))

#' @export
setAs("numeric", "lgCMatrix", function(from) as.csc.matrix(from, logical=TRUE))
#' @export
setAs("integer", "lgCMatrix", function(from) as.csc.matrix(from, logical=TRUE))
#' @export
setAs("logical", "lgCMatrix", function(from) as.csc.matrix(from, logical=TRUE))

#' @export
setAs("numeric", "ngCMatrix", function(from) as.csc.matrix(from, binary=TRUE))
#' @export
setAs("integer", "ngCMatrix", function(from) as.csc.matrix(from, binary=TRUE))
#' @export
setAs("logical", "ngCMatrix", function(from) as.csc.matrix(from, binary=TRUE))



#' @export
setAs("sparseMatrix", "dgTMatrix", function(from) as.coo.matrix(from))
#' @export
setAs("sparseMatrix", "lgTMatrix", function(from) as.coo.matrix(from, logical=TRUE))
#' @export
setAs("sparseMatrix", "ngTMatrix", function(from) as.coo.matrix(from, binary=TRUE))

#' @export
setAs("matrix", "dgTMatrix", function(from) as.coo.matrix(from))
#' @export
setAs("matrix", "lgTMatrix", function(from) as.coo.matrix(from, logical=TRUE))
#' @export
setAs("matrix", "ngTMatrix", function(from) as.coo.matrix(from, binary=TRUE))

#' @export
setAs("sparseVector", "dgTMatrix", function(from) as.coo.matrix(from))
#' @export
setAs("sparseVector", "lgTMatrix", function(from) as.coo.matrix(from, logical=TRUE))
#' @export
setAs("sparseVector", "ngTMatrix", function(from) as.coo.matrix(from, binary=TRUE))

#' @export
setAs("numeric", "dgTMatrix", function(from) as.coo.matrix(from))
#' @export
setAs("integer", "dgTMatrix", function(from) as.coo.matrix(from))
#' @export
setAs("logical", "dgTMatrix", function(from) as.coo.matrix(from))

#' @export
setAs("numeric", "lgTMatrix", function(from) as.coo.matrix(from, logical=TRUE))
#' @export
setAs("integer", "lgTMatrix", function(from) as.coo.matrix(from, logical=TRUE))
#' @export
setAs("logical", "lgTMatrix", function(from) as.coo.matrix(from, logical=TRUE))

#' @export
setAs("numeric", "ngTMatrix", function(from) as.coo.matrix(from, binary=TRUE))
#' @export
setAs("integer", "ngTMatrix", function(from) as.coo.matrix(from, binary=TRUE))
#' @export
setAs("logical", "ngTMatrix", function(from) as.coo.matrix(from, binary=TRUE))
