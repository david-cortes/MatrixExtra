
get_indices_integer <- function(i, max_i, index_names) {

    if (inherits(i, c("lsparseVector", "nsparseVector"))) {
        i <- sort_sparse_indices(i, copy=!getOption("MatrixExtra.inplace_sort", default=FALSE))
    }

    if (inherits(i, "lsparseVector")) {
        if (anyNA(i@x))
            stop("Subsetting with NA indices not supported.")
        i <- as(i, "nsparseVector")
    }

    if (inherits(i, "nsparseVector")) {
        if (i@length > max_i) {
            stop("Dimension of indexing vector is larger than matrix to subset.")
        } else if (i@length == max_i) {
            i <- i@i
        } else { ### mimic of R base's recycling
            full_repeats <- max_i %/% i@length
            remainder <- max_i - i@length*full_repeats
            if (remainder != 0)
                i <- repeat_indices_n_times(i@i, i[seq(1L, remainder)]@i, i@length, max_i)
            else
                i <- repeat_indices_n_times(i@i, integer(), i@length, max_i)
        }
        if (typeof(i) != "integer")
            i <- as.integer(i)
    }

    if (is.character(i)) {
        if (!is.null(index_names))
            i <- match(i, index_names)
        else
            i <- as.integer(i)
    }
    if (is.logical(i)) {
        if (length(i) != max_i) {
            i <- seq(1L, max_i)[i]
        } else {
            i <- which(i)
        }
    }
    ### TODO: maybe allow slicing with NAs
    if (anyNA(i))
        stop("Cannot slice matrix with NA indices.")
    if (typeof(i) != "integer")
        i <- as.integer(i)
    if (any(i <= 0))
        i <- seq(1L, max_i)[i]
    if(NROW(i) && any(i > max_i, na.rm=TRUE))
        stop("some of row subset indices are not present in matrix")
    return(as.integer(i))
}

ij_is_in_lower_triangle <- function(i, j) {
    return(i > j)
}

get_ij_properties <- function(x, i, j) {
    row_names <- rownames(x)
    col_names <- colnames(x)

    all_i <- FALSE
    all_j <- FALSE
    i_is_seq <- FALSE
    j_is_seq <- FALSE
    i_is_rev_seq <- FALSE
    j_is_rev_seq <- FALSE
    if (missing(j)) {
        all_j <- TRUE
        j <- seq_len(ncol(x))
        n_col <- ncol(x)
    } else {
        j <- get_indices_integer(j, ncol(x), col_names)
        if (length(j) == ncol(x) && j[1L] == 1L && j[length(j)] == ncol(x)) {
            if (check_is_seq(j))
                all_j <- TRUE
        } else {
            j_is_seq <- check_is_seq(j)
        }
        n_col <- length(j)
    }
    if (missing(i)) {
        i <- seq_len(nrow(x))
        all_i <- TRUE
        i_is_seq <- TRUE
        n_row <- nrow(x)
    } else {
        i <- get_indices_integer(i, nrow(x), row_names)
        i_is_seq <- check_is_seq(i)
        n_row <- length(i)

        if (i_is_seq && length(i) == nrow(x) && i[1L] == 1L && i[length(i)] == nrow(x)) {
            if (check_is_seq(i))
                all_i <- TRUE
        }
    }

    if (!all_i && !i_is_seq)
        i_is_rev_seq <- check_is_rev_seq(i)
    if (!all_j && !j_is_seq)
        j_is_rev_seq <- check_is_rev_seq(j)

    return(list(
        i = i, all_i = all_i, i_is_seq = i_is_seq, i_is_rev_seq = i_is_rev_seq,
        j = j, all_j = all_j, j_is_seq = j_is_seq, j_is_rev_seq = j_is_rev_seq,
        row_names = row_names, col_names = col_names,
        n_row = n_row, n_col = n_col
    ))
}

drop_slice <- function(x, drop) {
    
    if ((missing(drop) || isTRUE(drop)) && inherits(x, c("sparseMatrix", "float32"))) {
        if (nrow(x) == 1L || ncol(x) == 1L) {
            if (!(nrow(x) == 1L && ncol(x) == 1L) && !inherits(x, "float32") &&
                getOption("MatrixExtra.drop_sparse", default=FALSE)
            ) {
                x <- as(x, "sparseVector")
            } else {
                x <- as.vector(x)
            }
        }
    }

    if (inherits(x, "sparseMatrix")) {
        Dimnames <- list(NULL, NULL)
        if (NROW(x@Dimnames[[1L]]))
            Dimnames[[1L]] <- x@Dimnames[[1L]]
        if (NROW(x@Dimnames[[2L]]))
            Dimnames[[2L]] <- x@Dimnames[[2L]]
        x@Dimnames <- Dimnames
    }

    return(x)
}

#' @title Sparse Matrices Slicing
#' @description Natively slice CSR/CSC/COO matrices without changing the storage order.
#' @details \bold{Important:} When slicing sparse matrices with `drop=TRUE` (the default),
#' `Matrix` will drop 1-d matrices to \bold{dense} dense vectors, whereas
#' this package allows dropping them to either dense or \bold{sparse} vectors,
#' the latter of which is more efficient and is the default option.
#' 
#' The `drop` behavior can be changed back to dense vectors like `Matrix` does,
#' through \link{restore_old_matrix_behavior} or through the package options
#' (e.g. `options("MatrixExtra.drop_sparse" = FALSE)` - see \link{MatrixExtra-options}).
#' 
#' This package will override the subsetting methods from `Matrix` for all
#' sparse matrix types. It is usually much faster for all three storage orders (especially
#' CSR) but in some situations could end up being slightly slower. Be aware that, in the
#' case of COO matrices (a.k.a. "TsparseMatrix"), the resulting object will \bold{not}
#' have sorted indices, which `Matrix` will oftentimes do in addition to subsetting,
#' at a large speed penalty.
#' 
#' In general, it's much faster to select rows when the input is a CSR matrix ("RsparseMatrix"),
#' and much faster to select columns when the input is a CSC matrix ("CsparseMatrix").
#' Slicing COO matrices is typically not efficient, but could end up being faster when
#' the slice involves random rows and random columns with repeated entries.
#' @param x A sparse matrix to subset, in any format.
#' @param i row indices to subset.
#' @param j column indices to subset.
#' @param drop whether to simplify 1d matrix to a vector
#' @return A sparse matrix with the same storage order and dtype as `x`.
#' @name slice
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' m <- rsparsematrix(20, 20, 0.1, repr="R")
#' inherits(m[1:2, ], "RsparseMatrix")
#' inherits(m[1:2, 3:4], "RsparseMatrix")
#' inherits(as.coo.matrix(m)[1:2, 3:4], "TsparseMatrix")
#' inherits(as.csc.matrix(m)[1:2, 3:4], "CsparseMatrix")
#' 
#' ### New: slice with a sparse vector
#' m[as(c(TRUE,FALSE), "sparseVector"), ]
#' 
#' ### Important!!!
#' ### This differs from Matrix
#' set_new_matrix_behavior()
#' inherits(m[1,,drop=TRUE], "sparseVector")
#' 
#' ### To bring back the old behavior:
#' restore_old_matrix_behavior()
#' inherits(m[1,,drop=TRUE], "numeric")
NULL

subset_csr <- function(x, i, j, drop) {

    check_valid_matrix(x)

    if (missing(i) && missing(j)) {
        x <- drop_slice(x, drop)
        return(x)
    }

    if (!missing(i) && !missing(j) &&
        NROW(i) == 1L && NROW(j) == 1L &&
        typeof(i) %in% c("integer", "numeric") &&
        typeof(j) %in% c("integer", "numeric")
    ) {
        i <- as.integer(i)
        j <- as.integer(j)
        if (is.na(i) || is.na(j)) {
            res <- NA_real_
        } else if (i > nrow(x) || j > ncol(x)) {
            stop("Subscript out of bounds.")
        } else {

            ### TODO: add tests for these
            if (inherits(x, "triangularMatrix") &&
                ((x@uplo == "U" && ij_is_in_lower_triangle(i, j)) ||
                 (x@uplo == "L" && ij_is_in_lower_triangle(j, i)))
            ) {
                res <- 0
            } else if (.hasSlot(x, "diag") && x@diag == "U" && i == j) {
                res <- 1.
            } else {

                if (inherits(x, "symmetricMatrix")) {
                    if ((x@uplo == "U" && ij_is_in_lower_triangle(i, j)) ||
                        (x@uplo == "L" && ij_is_in_lower_triangle(j, i))
                    ) {
                        temp <- i
                        i <- j
                        j <- temp
                    }
                }

                if (inherits(x, "dsparseMatrix")) {
                    res <- extract_single_val_csr_numeric(x@p, x@j, x@x, i-1L, j-1L)
                } else if (inherits(x, "lsparseMatrix")) {
                    res <- as.logical(extract_single_val_csr_logical(x@p, x@j, x@x, i-1L, j-1L))
                } else if (inherits(x, "nsparseMatrix")) {
                    res <- extract_single_val_csr_binary(x@p, x@j, i-1L, j-1L)
                } else {
                    throw_internal_error()
                }
            }
        }
        if (missing(drop) || isTRUE(drop)) {
            if (!is.null(colnames(x)))
                names(res) <- colnames(x)[j]
            return(res)
        }

        if (inherits(x, "dsparseMatrix")) {
            out <- new("dgRMatrix")
        } else if (inherits(x, "lsparseMatrix")) {
            out <- new("lgRMatrix")
        } else if (inherits(x, "nsparseMatrix")) {
            out <- new("ngRMatrix")
        } else {
            throw_internal_error()
        }
        out@Dim <- c(1L, 1L)
        if (!is.null(x@Dimnames[[1L]]))
            out@Dimnames[[1L]] <- x@Dimnames[[1L]][i]
        if (!is.null(x@Dimnames[[2L]]))
            out@Dimnames[[2L]] <- x@Dimnames[[2L]][j]
        
        if (res == 0) {
            out@p <- c(0L, 0L)
            return(out)
        } else {
            out@p <- c(0L, 1L)
            out@j <- 0L
            if (.hasSlot(out, "x")) {
                if (inherits(out, "dsparseMatrix"))
                    out@x <- as.numeric(res)
                else
                    out@x <- as.logical(res)
            }
            return(out)
        }
    }

    ### TODO: should also add paths for sequences with a gap (e.g. 3, 5, 7, 9...)
    temp <- get_ij_properties(x, i, j)
    i <- temp$i
    j <- temp$j
    all_i <- temp$all_i
    all_j <- temp$all_j
    i_is_seq <- temp$i_is_seq
    j_is_seq <- temp$j_is_seq
    i_is_rev_seq <- temp$i_is_rev_seq
    j_is_rev_seq <- temp$j_is_rev_seq
    n_row <- temp$n_row
    n_col <- temp$n_col
    row_names <- temp$row_names
    col_names <- temp$col_names

    if (!NROW(i) || !NROW(j) || length(x@j) == 0L) {
        if (inherits(x, "dsparseMatrix")) {
            res <- new("dgRMatrix")
        } else if (inherits(x, "lsparseMatrix")) {
            res <- new("lgRMatrix")
        } else if (inherits(x, "nsparseMatrix")) {
            res <- new("ngRMatrix")
        } else {
            throw_internal_error()
        }
        res@p <- integer(NROW(i) + 1L)
        res@Dim <- c(NROW(i), NROW(j))

        row_names <- if(is.null(row_names) || !NROW(row_names)) NULL else row_names[i]
        col_names <- if(is.null(col_names) || !NROW(col_names)) NULL else col_names[j]
        res@Dimnames <- list(row_names, col_names)
        res <- drop_slice(res, drop)
        return(res)
    }

    if (all_i && all_j) {
        return(drop_slice(x, drop))
    }

    if (!all_i && !i_is_seq)
        i_is_rev_seq <- check_is_rev_seq(i)
    if (!all_j && !j_is_seq)
        j_is_rev_seq <- check_is_rev_seq(j)

    if (i_is_rev_seq && j_is_rev_seq && length(i) == nrow(x) && length(j) == ncol(x)) {
        if (inherits(x, "dsparseMatrix")) {
            temp <- reverse_rows_numeric(x@p, x@j, x@x)
            reverse_columns_inplace_numeric(temp$indptr, temp$indices, temp$values, ncol(x))
        } else if (inherits(x, "lsparseMatrix")) {
            temp <- reverse_rows_logical(x@p, x@j, x@x)
            reverse_columns_inplace_logical(temp$indptr, temp$indices, temp$values, ncol(x))
        } else if (inherits(x, "nsparseMatrix")) {
            temp <- reverse_rows_binary(x@p, x@j)
            reverse_columns_inplace_binary(temp$indptr, temp$indices, ncol(x))
        } else {
            throw_internal_error()
        }

        X_attr <- attributes(x)
        X_attr$p <- temp$indptr
        X_attr$j <- temp$indices
        if ("x" %in% names(X_attr))
            X_attr$x <- temp$values
        if (!is.null(X_attr$Dimnames[[1L]]))
            X_attr$Dimnames[[1L]] <- rev(X_attr$Dimnames[[1L]])
        if (!is.null(X_attr$Dimnames[[2L]]))
            X_attr$Dimnames[[2L]] <- rev(X_attr$Dimnames[[2L]])
        if ("uplo" %in% names(X_attr))
            X_attr$uplo = ifelse(X_attr$uplo == "U", "L", "U")
        attributes(x) <- X_attr
        return(drop_slice(x, drop))
    }

    if (inherits(x, c("symmetricMatrix", "triangularMatrix")) && !(length(x@j) == 0L))
        x <- as.csr.matrix(x, logical=inherits(x, "lsparseMatrix"), binary=inherits(x, "nsparseMatrix"))
    has_x <- .hasSlot(x, "x")

    if (i_is_seq && all_j) {
        first <- x@p[i[1L]] + 1L
        last <- x@p[i[n_row] + 1L]
        indptr <- x@p[seq(i[1L], i[n_row]+1L)] - x@p[i[1L]]
        col_indices <- x@j[first:last]
        if (has_x)
            x_values <- x@x[first:last]
    } else if (i_is_rev_seq && all_j) {
        first <- x@p[i[n_row]] + 1L
        last <- x@p[i[1L] + 1L]
        indptr <- x@p[seq(i[n_row], i[1L]+1L)] - x@p[i[n_row]]
        col_indices <- x@j[first:last]
        if (has_x)
            x_values <- x@x[first:last]
        if (!has_x) {
            temp <- reverse_rows_binary(indptr, col_indices)
        } else if (typeof(x_values) == "logical") {
            temp <- reverse_rows_logical(indptr, col_indices, x_values)
        } else {
            temp <- reverse_rows_numeric(indptr, col_indices, as.numeric(x_values))
        }
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    } else if (!i_is_seq && all_j) {
        if (inherits(x, "dsparseMatrix")) {
            temp <- copy_csr_rows_numeric(x@p, x@j, x@x, i-1L)
        } else if (inherits(x, "lsparseMatrix")) {
            temp <- copy_csr_rows_logical(x@p, x@j, x@x, i-1L)
        } else if (inherits(x, "nsparseMatrix")) {
            temp <- copy_csr_rows_binary(x@p, x@j, i-1L)
        } else {
            throw_internal_error()
        }
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    } else if (j_is_seq) {
        if (inherits(x, "dsparseMatrix")) {
            temp <- copy_csr_rows_col_seq_numeric(x@p, x@j, x@x, i-1L, j, TRUE)
        } else if (inherits(x, "lsparseMatrix")) {
            temp <- copy_csr_rows_col_seq_logical(x@p, x@j, x@x, i-1L, j, TRUE)
        } else if (inherits(x, "nsparseMatrix")) {
            temp <- copy_csr_rows_col_seq_binary(x@p, x@j, i-1L, j, TRUE)
        } else {
            throw_internal_error()
        }
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    } else if (j_is_rev_seq) {
        new_ncol <- j[1L] - j[length(j)] + 1L
        if (inherits(x, "dsparseMatrix")) {
            temp <- copy_csr_rows_col_seq_numeric(x@p, x@j, x@x, i-1L, j, TRUE)
            reverse_columns_inplace_numeric(temp$indptr, temp$indices, temp$values, new_ncol)
        } else if (inherits(x, "lsparseMatrix")) {
            temp <- copy_csr_rows_col_seq_logical(x@p, x@j, x@x, i-1L, j, TRUE)
            reverse_columns_inplace_logical(temp$indptr, temp$indices, temp$values, new_ncol)
        } else if (inherits(x, "nsparseMatrix")) {
            temp <- copy_csr_rows_col_seq_binary(x@p, x@j, i-1L, j, TRUE)
            reverse_columns_inplace_binary(temp$indptr, temp$indices, new_ncol)
        } else {
            throw_internal_error()
        }
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    } else {
        ### FAIL: this ended up being much slower than the original Matrix
        ### when the columns don't get pre-discarded through the
        ### condition min(j) <= col < max(j)
        if (inherits(x, "dsparseMatrix")) {
            temp <- copy_csr_arbitrary_numeric(x@p, x@j, x@x, i-1L, j-1L)
        } else if (inherits(x, "lsparseMatrix")) {
            temp <- copy_csr_arbitrary_logical(x@p, x@j, x@x, i-1L, j-1L)
        } else if (inherits(x, "nsparseMatrix")) {
            temp <- copy_csr_arbitrary_binary(x@p, x@j, i-1L, j-1L)
        } else {
            throw_internal_error()
        }
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    }

    res <- new(class(x)[1L])
    res@p <- indptr
    res@j <- col_indices
    if (has_x) {
        if (inherits(x, "lsparseMatrix")) {
            res@x <- as.logical(x_values)
        } else {
            res@x <- x_values
        }
    }
    res@Dim <- c(n_row, n_col)

    row_names <- if (is.null(row_names) || !NROW(row_names)) NULL else row_names[i]
    col_names <- if (is.null(col_names) || !NCOL(col_names)) NULL else col_names[j]
    res@Dimnames <- list(row_names, col_names)

    res <- drop_slice(res, drop)
    return(res)
}

#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="index", j="index", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="missing", j="index", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="index", j="missing", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="missing", j="missing", drop="logical"), subset_csr)

#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="index", j="index", drop="missing"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="missing", j="index", drop="missing"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="index", j="missing", drop="missing"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="RsparseMatrix", i="missing", j="missing", drop="missing"), subset_csr)

subset_generic_with_vector <- function(x, i, j, drop) {
    if (inherits(x, "CsparseMatrix"))
        return(t_shallow(subset_csr(t_shallow(x), j, i, drop)))
    if (inherits(x, "RsparseMatrix"))
        return(subset_csr(x, i, j, drop))

    if (missing(i) && missing(j))
        return(drop_slice(x, drop))
    if (!missing(i))
        i <- get_indices_integer(i, NROW(x), rownames(x))
    if (!missing(j))
        j <- get_indices_integer(j, NCOL(x), colnames(x))

    if (!is.null(dim(x))) {
        if (!missing(i) && !missing(j)) {
            return(x[i, j, drop=drop])
        } else if (!missing(i)) {
            return(x[i, , drop=drop])
        } else {
            return(x[, j, drop=drop])
        }
    } else {
        if (missing(i))
            return(x)
        return(x[i])
    }
}

#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="nsparseVector", j="nsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="missing", j="nsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="nsparseVector", j="missing", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="index", j="nsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="nsparseVector", j="index", drop="logical"), subset_generic_with_vector)


#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="nsparseVector", j="nsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="missing", j="nsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="nsparseVector", j="missing", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="index", j="nsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="nsparseVector", j="index", drop="missing"), subset_generic_with_vector)


#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="lsparseVector", j="lsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="missing", j="lsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="lsparseVector", j="missing", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="index", j="lsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="lsparseVector", j="index", drop="logical"), subset_generic_with_vector)


#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="lsparseVector", j="lsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="missing", j="lsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="lsparseVector", j="missing", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="index", j="lsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="ANY", i="lsparseVector", j="index", drop="missing"), subset_generic_with_vector)


### The method  here is usually faster than CSC slicing from Matrix, so it will replace it too

subset_csc_masked <- function(x, i, j, drop) {
    x <- subset_csr(t_shallow(x), j, i, drop)
    if (inherits(x, "sparseMatrix"))
        x <- t_shallow(x)
    return(x)
}

#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="index", j="index", drop="logical"), subset_csc_masked)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="missing", j="index", drop="logical"), subset_csc_masked)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="index", j="missing", drop="logical"), subset_csc_masked)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="missing", j="missing", drop="logical"), subset_csc_masked)

#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="index", j="index", drop="missing"), subset_csc_masked)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="missing", j="index", drop="missing"), subset_csc_masked)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="index", j="missing", drop="missing"), subset_csc_masked)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="CsparseMatrix", i="missing", j="missing", drop="missing"), subset_csc_masked)
