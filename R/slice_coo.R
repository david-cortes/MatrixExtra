### Note: this method is generally not any faster than the one in Matrix,
### but having it redone here allows keeping the drop-to-sparse option
### for all input types.
### Additionally, there are issues when trying to use the function that Matrix
### will call to slice COO matrices due to pontential problems with circular
### references after having masked methods for CSC too.

subset_coo <- function(x, i, j, drop) {
    
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

                ix_take <- slice_coo_single(x@i, x@j, i-1L, j-1L)
                if (ix_take == 0) {
                    res <- 0.
                } else {
                    if (.hasSlot(x, "x"))
                        res <- x@x[ix_take]
                    else
                        res <- 1.
                }
            }
        }
        if (missing(drop) || isTRUE(drop)) {
            if (!is.null(colnames(x)))
                names(res) <- colnames(x)[j]
            return(res)
        }

        if (inherits(x, "dsparseMatrix")) {
            out <- new("dgTMatrix")
        } else if (inherits(x, "lsparseMatrix")) {
            out <- new("lgTMatrix")
        } else if (inherits(x, "nsparseMatrix")) {
            out <- new("ngTMatrix")
        } else {
            stop("Internal error. Please open an issue in GitHub.")
        }
        out@Dim <- c(1L, 1L)
        if (!is.null(x@Dimnames[[1L]]))
            out@Dimnames[[1L]] <- x@Dimnames[[1L]][i]
        if (!is.null(x@Dimnames[[2L]]))
            out@Dimnames[[2L]] <- x@Dimnames[[2L]][j]
        
        if (is.na(res) || res != 0) {
            out@i <- 0L
            out@j <- 0L
            if (inherits(out, "dsparseMatrix"))
                out@x <- as.numeric(res)
            else
                out@x <- as.logical(res)
        }
        return(out)
    }


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
    

    if (!NROW(i) || !NROW(j) || length(x@i) == 0L) {
        if (inherits(x, "dsparseMatrix")) {
            res <- new("dgTMatrix")
        } else if (inherits(x, "lsparseMatrix")) {
            res <- new("lgTMatrix")
        } else if (inherits(x, "nsparseMatrix")) {
            res <- new("ngTMatrix")
        } else {
            stop("Internal error. Please open an issue in GitHub.")
        }
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

    if (i_is_rev_seq && j_is_rev_seq && length(i) == nrow(x) && length(j) == ncol(x)) {
        X_attr <- attributes(x)
        X_attr$i <- (nrow(x) - 1L) - X_attr$i
        X_attr$j <- (ncol(x) - 1L) - X_attr$j
        if (!is.null(X_attr$Dimnames[[1L]]))
            X_attr$Dimnames[[1L]] <- rev(X_attr$Dimnames[[1L]])
        if (!is.null(X_attr$Dimnames[[2L]]))
            X_attr$Dimnames[[2L]] <- rev(X_attr$Dimnames[[2L]])
        if ("uplo" %in% names(X_attr))
            X_attr$uplo = ifelse(X_attr$uplo == "U", "L", "U")
        attributes(x) <- X_attr
        return(drop_slice(x, drop))
    }

    if (inherits(x, c("symmetricMatrix", "triangularMatrix")))
        x <- as.coo.matrix(x, logical=inherits(x, "lsparseMatrix"), binary=inherits(x, "nsparseMatrix"))
    has_x <- .hasSlot(x, "x")

    if (inherits(x, "dsparseMatrix")) {
        temp <- slice_coo_arbitrary_numeric(
            x@i, x@j, x@x,
            i, j,
            all_i, all_j,
            i_is_seq, j_is_seq,
            i_is_rev_seq, j_is_rev_seq,
            nrow(x), ncol(x)
        )
    } else if (inherits(x, "lsparseMatrix")) {
        temp <- slice_coo_arbitrary_logical(
            x@i, x@j, x@x,
            i, j,
            all_i, all_j,
            i_is_seq, j_is_seq,
            i_is_rev_seq, j_is_rev_seq,
            nrow(x), ncol(x)
        )
    } else if (inherits(x, "nsparseMatrix")) {
        temp <- slice_coo_arbitrary_binary(
            x@i, x@j,
            i, j,
            all_i, all_j,
            i_is_seq, j_is_seq,
            i_is_rev_seq, j_is_rev_seq,
            nrow(x), ncol(x)
        )
    } else {
        throw_internal_error()
    }

    res <- new(class(x)[1L])
    res@i <- temp$ii
    res@j <- temp$jj
    if (.hasSlot(res, "x"))
        res@x <- temp$xx

    res@Dim <- c(n_row, n_col)

    row_names <- if (is.null(row_names) || !NROW(row_names)) NULL else row_names[i]
    col_names <- if (is.null(col_names) || !NCOL(col_names)) NULL else col_names[j]
    res@Dimnames <- list(row_names, col_names)

    res <- drop_slice(res, drop)
    return(res)
} 


#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="index", j="index", drop="logical"), subset_coo)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="missing", j="index", drop="logical"), subset_coo)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="index", j="missing", drop="logical"), subset_coo)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="missing", j="missing", drop="logical"), subset_coo)

#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="index", j="index", drop="missing"), subset_coo)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="missing", j="index", drop="missing"), subset_coo)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="index", j="missing", drop="missing"), subset_coo)
#' @rdname slice
#' @export
setMethod(`[`, signature(x="TsparseMatrix", i="missing", j="missing", drop="missing"), subset_coo)
