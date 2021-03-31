#' @title Concatenate inputs by rows into a CSR matrix
#' @description Concatenate two or more matrices and/or vectors by rows, giving a CSR matrix
#' as result.
#'
#' This is aimed at concatenating several CSR matrices or sparse vectors at a time,
#' as it will be faster than calling `rbind` which will only concatenate one at a
#' time, resulting in unnecessary allocations.
#' @param ... Inputs to concatenate. The function is aimed at CSR matrices (`dgRMatrix`,
#' `ngRMatrix`, `lgRMatrix`) and sparse vectors (`sparseVector`). It will work with other classes
#' (such as `dgCMatrix`) but will not be as efficient.
#' @returns A CSR matrix (class `dgRMatrix`, `lgRMatrix`, or `ngRMatrix` depending on the inputs) with
#' the inputs concatenated by rows.
#' @details This function will not preserve the column names, if any were present.
#' @seealso \link{rbind2-method}
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' v <- as(1:10, "sparseVector")
#' rbind_csr(v, v, v)
#'
#' X <- matrix(1:20, nrow=2)
#' rbind_csr(X, v)
#' @export
rbind_csr <- function(...) {

    binary_types <- c("nsparseMatrix", "nsparseVector")
    logical_types <- c("lsparseMatrix", "lsparseVector", "logical")
    cast_csr_same <- function(x) as.csr.matrix(x, binary=inherits(x, binary_types), logical=inherits(x, logical_types))
    cast_if_not_csr <- function(x) {
        if (inherits(x, c("dgRMatrix", "ngRMatrix", "lgRMatrix", "sparseVector"))) {
            return(x)
        } else {
            return(cast_csr_same(x))
        }
    }
    args <- lapply(list(...), cast_if_not_csr)

    if (length(args) == 0L) {
        return(new("dgRMatrix"))
    } else if (length(args) == 1L) {
        return(cast_csr_same(args[[1L]]))
    } else if (length(args) == 2L) {
        return(rbind2(args[[1L]], args[[2L]]))
    }

    is_binary <- sapply(args, function(x) inherits(x, binary_types))
    is_logical <- sapply(args, function(x) inherits(x, logical_types))
    is_numeric <- !is_binary & !is_logical
    if (any(is_numeric)) {
        out <- new("dgRMatrix")
    } else if (any(is_logical)) {
        out <- new("lgRMatrix")
    } else {
        out <- new("ngRMatrix")
    }

    nrows <- sum(sapply(args, function(x) ifelse(inherits(x, "sparseMatrix"), nrow(x), 1L)))
    ncols <- max(sapply(args, function(x) ifelse(inherits(x, "sparseMatrix"), ncol(x), length(x))))
    nnz <- sum(sapply(args, function(x) ifelse(inherits(x, "sparseMatrix"), length(x@j), length(x@i))))
    if (nrows >= .Machine$integer.max)
        stop("Result has too many rows for R to handle.")
    if (nnz >= .Machine$integer.max)
        stop("Result has too many non-zero entries for R to handle.")

    out@p <- integer(nrows + 1L)
    out@j <- integer(nnz)
    if (inherits(out, "dgRMatrix")) {
        out@x <- numeric(nnz)
    } else if (inherits(out, "lgRMatrix")) {
        out@x <- logical(nnz)
    }
    out@Dim <- as.integer(c(nrows, ncols))
    out@Dimnames <- list(NULL, NULL)
    if (!nrows || !ncols)
        return(out)

    out <- concat_csr_batch(args, out)

    if (any(sapply(args, function(x) !is.null(rownames(x))))) {
        rownames(out) <- Reduce(c, lapply(args, function(x) {
            if (!is.null(rownames(x))) {
                return(as.character(rownames(x)))
            } else {
                if (inherits(x, "sparseVector")) {
                    return("")
                } else {
                    return(rep("", nrow(x)))
                }
            }
        }))
    }

    return(out)
}

#' @name rbind2-method
#' @title Concatenate CSR matrices by rows
#' @description `rbind2` method for the sparse matrix and sparse vector classes from `Matrix`,
#' taking the most efficient route for the concatenation according to the input types.
#' @param x First matrix to concatenate.
#' @param y Second matrix to concatenate.
#' @return Either a CSR matrix or a COO matrix depending on the input types.
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(3, 4, .3)
#' X <- as(X, "RsparseMatrix")
#' inherits(rbind2(X, X), "dgRMatrix")
#' inherits(rbind(X, X, as.csc.matrix(X), X), "dgRMatrix")
#' inherits(rbind2(as.coo.matrix(X), as.coo.matrix(X)), "dgTMatrix")
NULL

concat_dimname <- function(nm1, nm2, nrow1, nrow2) {
    if (is.null(nm1) && is.null(nm2)) {
        return(NULL)
    } else if (!is.null(nm1) && !is.null(nm2)) {
        return(c(nm1, nm2))
    }

    if (is.null(nm1)) {
        nm1 <- as.character(seq(1, nrow1))
    }
    if (is.null(nm2)) {
        nm2 <- as.character(seq(nrow1 + 1L, nrow1 + nrow2))
    }
    return(concat_dimname(nm1, nm2, nrow1, nrow2))
}

concat_dimnames <- function(mat1, mat2) {
    dim1 <- concat_dimname(mat1@Dimnames[[1L]], mat2@Dimnames[[1L]], nrow(mat1), nrow(mat2))

    if (!is.null(mat1@Dimnames[[2L]])) {
        dim2 <- mat1@Dimnames[[2L]]
        if (mat1@Dim[2L] < mat2@Dim[2L]) {
            if (!is.null(mat2@Dimnames[[2L]])) {
                dim2 <- c(dim2, mat2@Dimnames[[2L]][seq(length(dim2)+1,length(mat2@Dimnames[[2L]]))])
            } else {
                dim2 <- c(dim2, rep("", length(mat2@Dimnames[[2L]]) - length(dim2)))
            }
        }
    } else if (!is.null(mat2@Dimnames[[2L]])) {
        dim2 <- mat2@Dimnames[[2L]]
    } else {
        dim2 <- NULL
    }

    if (!is.null(dim1) && typeof(dim1) != "character") {
        dim1 <- as.character(dim1)
        largest <- max(mat1@Dim[1L], mat2@Dim[1L])
        if (length(dim1) < largest)
            dim1 <- c(dim1, rep("", largest - length(dim1)))
    }
    if (!is.null(dim2) && typeof(dim2) != "character") {
        dim2 <- as.character(dim2)
        largest <- max(mat1@Dim[2L], mat2@Dim[2L])
        if (length(dim2) < largest)
            dim2 <- c(dim2, rep("", largest - length(dim2)))
    }

    return(list(dim1, dim2))
}

concat_as_numeric <- function(v1, v2) {
    if (typeof(v1) != "double") mode(v1) <- "double"
    if (typeof(v2) != "double") mode(v2) <- "double"
    return(c(v1, v2))
}

concat_as_logical <- function(v1, v2) {
    if (typeof(v1) != "logical") mode(v1) <- "logical"
    if (typeof(v2) != "logical") mode(v2) <- "logical"
    return(c(v1, v2))
}

rbind2_csr <- function(x, y, out) {
    Dim <- c(x@Dim[1L] + y@Dim[1L], max(x@Dim[2L], y@Dim[2L]))
    if (Dim[2L] >= .Machine$integer.max)
        stop("Resulting matrix has too many rows for R to handle.")
    out@Dim <- as.integer(Dim)
    out@Dimnames <- concat_dimnames(x, y)
    out@p <- concat_indptr2(x@p, y@p)
    out@j <- c(x@j, y@j)
    return(out)
}

rbind2_dgr <- function(x, y) {
    out <- new("dgRMatrix")
    out <- rbind2_csr(x, y, out)
    if (.hasSlot(x, "x") && .hasSlot(y, "x")) {
        out@x <- concat_as_numeric(x@x, y@x)
    } else if (.hasSlot(x, "x")) {
        out@x <- concat_as_numeric(x@x, rep(1., length(y@j)))
    } else if (.hasSlot(y, "x")) {
        out@x <- concat_as_numeric(rep(1., length(x@j)), y@x)
    } else {
        out@x <- rep(1., length(x@j) + length(y@j))
    }
    return(out)
}

rbind2_lgr <- function(x, y) {
    out <- new("lgRMatrix")
    out <- rbind2_csr(x, y, out)
    if (.hasSlot(x, "x") && .hasSlot(y, "x")) {
        out@x <- concat_as_logical(x@x, y@x)
    } else if (.hasSlot(x, "x")) {
        out@x <- concat_as_logical(x@x, rep(TRUE, length(y@j)))
    } else if (.hasSlot(y, "x")) {
        out@x <- concat_as_logical(rep(TRUE, length(x@j)), y@x)
    } else {
        out@x <- rep(TRUE, length(x@j) + length(y@j))
    }
    return(out)
}

rbind2_ngr <- function(x, y) {
    out <- new("ngRMatrix")
    out <- rbind2_csr(x, y, out)
    return(out)
}

rbind2_generic <- function(x, y) {
    binary_types <- c("nsparseMatrix", "nsparseVector")
    logical_types <- c("lsparseMatrix", "lsparseVector")
    x_is_binary <- inherits(x, binary_types)
    x_is_logical <- inherits(x, logical_types)
    y_is_binary <- inherits(y, binary_types)
    y_is_logical <- inherits(y, logical_types)

    if (inherits(x, "symmetricMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.csr.matrix(x, logical=x_is_logical, binary=x_is_binary)
    if (inherits(x, "symmetricMatrix") || (.hasSlot(y, "diag") && y@diag != "N"))
        y <- as.csr.matrix(y, logical=y_is_logical, binary=y_is_binary)

    if (x_is_binary && y_is_binary) {
        return(rbind2_ngr(as.csr.matrix(x, binary=TRUE), as.csr.matrix(y, binary=TRUE)))
    } else if ((x_is_binary || x_is_logical) && (y_is_binary || y_is_logical)) {
        return(rbind2_lgr(as.csr.matrix(x, logical=TRUE), as.csr.matrix(y, logical=TRUE)))
    } else {
        return(rbind2_dgr(as.csr.matrix(x), as.csr.matrix(y)))
    }
}

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="RsparseMatrix", y="RsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="sparseVector", y="RsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="RsparseMatrix", y="sparseVector"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="sparseVector", y="sparseVector"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="CsparseMatrix", y="CsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="sparseVector", y="CsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="CsparseMatrix", y="sparseVector"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="CsparseMatrix", y="CsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="RsparseMatrix", y="CsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="CsparseMatrix", y="RsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="RsparseMatrix", y="TsparseMatrix"), rbind2_generic)

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="TsparseMatrix", y="RsparseMatrix"), rbind2_generic)

### TODO: add tests for the ones below

rbind2_tri <- function(x, y) {
    x_is_binary <- inherits(x, "nsparseMatrix")
    x_is_logical <- inherits(x, "lsparseMatrix")
    y_is_binary <- inherits(y, "nsparseMatrix")
    y_is_logical <- inherits(y, "lsparseMatrix")

    if (inherits(x, "symmetricMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x, logical=x_is_logical, binary=x_is_binary)
    if (inherits(x, "symmetricMatrix") || (.hasSlot(y, "diag") && y@diag != "N"))
        y <- as.coo.matrix(y, logical=y_is_logical, binary=y_is_binary)

    if (x_is_binary && y_is_binary) {
        out <- new("ngTMatrix")
    } else if ((x_is_binary || x_is_logical) && (y_is_binary || y_is_logical)) {
        out <- new("lgTMatrix")
    } else {
        out <- new("dgTMatrix")
    }

    out@Dim <- as.integer(c(nrow(x) + nrow(y), max(ncol(x), ncol(y))))
    out@Dimnames <- concat_dimnames(x, y)
    out@i <- c(x@i, y@i)
    out@j <- c(x@j, i@j)

    if (inherits(out, "dsparseMatrix")) {
        if (.hasSlot(x, "x") && .hasSlot(y, "x")) {
            out@x <- concat_as_numeric(x@x, y@x)
        } else if (.hasSlot(x, "x")) {
            out@x <- concat_as_numeric(x@x, rep(1, length(y@i)))
        } else {
            out@x <- concat_as_numeric(rep(1, length(x@i)), y@x)
        }
    } else if (inherits(out, "lsparseMatrix")) {
        if (.hasSlot(x, "x") && .hasSlot(y, "x")) {
            out@x <- concat_as_logical(x@x, y@x)
        } else if (.hasSlot(x, "x")) {
            out@x <- concat_as_logical(x@x, rep(TRUE, length(y@i)))
        } else {
            out@x <- concat_as_logical(rep(TRUE, length(x@i)), y@x)
        }
    }

    return(out)
}

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="TsparseMatrix", y="TsparseMatrix"), rbind2_tri)

rbind2_tri_vec <- function(x, v, x_is_first) {
    x_is_binary <- inherits(x, "nsparseMatrix")
    x_is_logical <- inherits(x, "lsparseMatrix")
    v_is_binary <- inherits(v, "nsparseVector")
    v_is_logical <- inherits(v, "lsparseVector")

    if (inherits(x, "symmetricMatrix") || (.hasSlot(x, "diag") && x@diag != "N"))
        x <- as.coo.matrix(x, logical=x_is_logical, binary=x_is_binary)

    if (x_is_binary && v_is_binary) {
        out <- new("ngTMatrix")
    } else if ((x_is_binary || x_is_logical) && (v_is_binary || v_is_logical)) {
        out <- new("lgTMatrix")
    } else {
        out <- new("dgTMatrix")
    }

    out@Dim <- as.integer(c(nrow(x) + 1L, max(ncol(x), v@length)))
    if (x_is_first)
        out@Dimnames <- concat_dimnames(x, v)
    else
        out@Dimnames <- concat_dimnames(v, x)
    out@i <- c(x@i, v@i - 1L)
    out@j <- c(x@j, i@j)

    if (inherits(out, "dsparseMatrix")) {
        if (.hasSlot(x, "x") && .hasSlot(v, "x")) {
            out@x <- concat_as_numeric(x@x, v@x)
        } else if (.hasSlot(x, "x")) {
            out@x <- concat_as_numeric(x@x, rep(1, length(v@i)))
        } else {
            out@x <- concat_as_numeric(rep(1, length(x@i)), v@x)
        }
    } else if (inherits(out, "lsparseMatrix")) {
        if (.hasSlot(x, "x") && .hasSlot(v, "x")) {
            out@x <- concat_as_logical(x@x, v@x)
        } else if (.hasSlot(x, "x")) {
            out@x <- concat_as_logical(x@x, rep(TRUE, length(v@i)))
        } else {
            out@x <- concat_as_logical(rep(TRUE, length(x@i)), v@x)
        }
    }

    return(out)
}

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="TsparseMatrix", y="sparseVector"), function(x, y) rbind2_tri_vec(x, y))

#' @rdname rbind2-method
#' @export
setMethod("rbind2", signature(x="sparseVector", y="TsparseMatrix"), function(x, y) rbind2_tri_vec(y, x, FALSE))
