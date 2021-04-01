
get_indices_integer <- function(i, max_i, index_names) {

    if (inherits(i, "nsparseVector")) {
        if (i@length > max_i) {
            stop("Dimension of indexing vector is larger than matrix to subset.")
        } else if (i@length == max_i) {
            i <- i@i
        } else { ### mimic of R base's recycling
            full_repeats <- max_i %/% i@length
            remainder <- max_i - i@length*full_repeats
            i <- repeat_indices_n_times(i@i, i[seq(1L, remainder)]@i, i@length, max_i)
        }
        if (typeof(i) != "integer")
            i <- as.integer(i)
    }

    if (is.numeric(i)) i <- as.integer(i)
    if (is.character(i)) i <- match(i, index_names)
    if (is.logical(i)) {
        if (length(i) != max_i) {
            i <- seq(1L, max_i)[i]
        } else {
            i <- which(i)
        }
    }
    if (any(i < 0))
        i <- seq(1L, max_i)[i]
    if(anyNA(i) || any(i >    max_i, na.rm=TRUE))
        stop("some of row subset indices are not present in matrix")
    return(as.integer(i))
}

get_x_values <- function(Mat) {
    if (.hasSlot(Mat, "x")) {
        return(as.numeric(Mat@x))
    } else {
        return(numeric())
    }
}

#' @title CSR Matrices Slicing
#' @description Natively slice CSR matrices without converting them to triplets/CSC.
#' @param x input `RsparseMatrix`
#' @param i row indices to subset
#' @param j column indices to subset
#' @param drop whether to simplify 1d matrix to a (dense) vector
#' @return
#' An `RsparseMatrix`.
#' @name slice
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' # dgCMatrix - CSC
#' m <- rsparsematrix(20, 20, 0.1)
#' # make CSR
#' m <- as(m, "RsparseMatrix")
#' inherits(m[1:2, ], "RsparseMatrix")
#' inherits(m[1:2, 3:4], "RsparseMatrix")
NULL

subset_csr <- function(x, i, j, drop=TRUE) {

    if (inherits(x, c("symmetricMatrix", "triangularMatrix")))
        x <- as.csr.matrix(x)
    has_x <- .hasSlot(x, "x")

    if (missing(i) && missing(j)) {
        if (missing(drop) || isTRUE(drop)) {
            if (nrow(x) == 1L || ncol(x) == 1L)
                x <- as.vector(x)
        }
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
            if (.hasSlot(x, "x"))
                res <- extract_single_val_csr_numeric(x@p, x@j, x@x, i-1L, j-1L)
            else
                res <- extract_single_val_csr_binary(x@p, x@j, i-1L, j-1L)
        }
        if (missing(drop) || isTRUE(drop)) {
            if (!is.null(colnames(x)))
                names(res) <- colnames(x)[j]
            return(res)
        }

        out <- new("dgRMatrix")
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
            out@x <- res
            return(out)
        }
    }

    row_names <- rownames(x)
    col_names <- colnames(x)

    all_i <- FALSE
    all_j <- FALSE
    i_is_seq <- FALSE
    j_is_seq <- FALSE
    if (missing(j)) {
        all_j <- TRUE
        j <- seq_len(ncol(x))
        n_col <- ncol(x)
    } else {
        j <- get_indices_integer(j, ncol(x), col_names)
        if (length(j) == ncol(x) && j[1L] == 1L && j[length(j)] == ncol(x)) {
            if (all(j == seq(1L, ncol(x))))
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
        # convert integer/numeric/logical/character indices to integer indices
        # also takes care of negative indices
        i <- get_indices_integer(i, nrow(x), row_names)
        i_is_seq <- check_is_seq(i)
        n_row <- length(i)

        if (all_j && i_is_seq && length(i) == nrow(x) && i[1L] == 1L && i[length(i)] == nrow(x)) {
            if (all(i == seq(1L, nrow(x))))
                all_i <- TRUE
        }
    }

    if (!NROW(i) || !NROW(j)) {
        res <- new(class(x))
        res@p <- integer(NROW(i) + 1L)
        res@Dim <- c(NROW(i), NROW(j))

        row_names <- if(is.null(row_names) || !NROW(row_names)) NULL else row_names[i]
        col_names <- if(is.null(col_names) || !NROW(col_names)) NULL else col_names[j]
        res@Dimnames <- list(row_names, col_names)
        return(res)
    }

    if (all_i && all_j) {
        return(x)
    } else if (length(x@j) == 0L) {
        indptr <- integer(n_row + 1)
        col_indices <- integer()
        x_values <- numeric()
    } else if (i_is_seq && all_j) {
        first <- x@p[i[1L]] + 1L
        last <- x@p[i[n_row] + 1L]
        indptr <- x@p[seq(i[1L], i[n_row]+1L)] - x@p[i[1L]]
        col_indices <- x@j[first:last]
        if (has_x)
            x_values <- x@x[first:last]
    } else if (!i_is_seq && all_j) {
        temp <- copy_csr_rows(x@p, x@j, get_x_values(x), i-1L)
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    } else if (j_is_seq) {
        temp <- copy_csr_rows_col_seq(x@p, x@j, get_x_values(x), i-1L, j-1L)
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    } else {
        temp <- copy_csr_arbitrary(x@p, x@j, get_x_values(x), i-1L, j-1L)
        indptr <- temp$indptr
        col_indices <- temp$indices
        if (has_x)
            x_values <- temp$values
    }

    res <- new(class(x)[1L])
    res@p <- indptr
    res@j <- col_indices
    if (has_x) {
        if (inherits(x, "lgRMatrix")) {
            res@x <- as.logical(x_values)
        } else {
            res@x <- x_values
        }
    }
    res@Dim <- c(n_row, n_col)

    row_names <- if (is.null(row_names) || !NROW(row_names)) NULL else row_names[i]
    col_names <- if (is.null(col_names) || !NCOL(col_names)) NULL else col_names[j]
    res@Dimnames <- list(row_names, col_names)

    if(isTRUE(drop) && (n_row == 1L || n_col == 1L))
        res <- as.vector(res)
    return(res)
}

#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "index", j <- "index", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "missing", j <- "index", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "index", j <- "missing", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "missing", j <- "missing", drop="logical"), subset_csr)

#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "nsparseVector", j <- "nsparseVector", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "missing", j <- "nsparseVector", drop="logical"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "nsparseVector", j <- "missing", drop="logical"), subset_csr)

#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "nsparseVector", j <- "nsparseVector", drop="missing"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "missing", j <- "nsparseVector", drop="missing"), subset_csr)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "RsparseMatrix", i <- "nsparseVector", j <- "missing", drop="missing"), subset_csr)

subset_generic_with_vector <- function(x, i, j, drop) {
    if (missing(i) && missing(j))
        return(x)
    if (!missing(i))
        i <- get_indices_integer(i, nrow(x), rownames(x))
    if (!missing(j))
        j <- get_indices_integer(j, ncol(x), colnames(x))

    if (!missing(i) && !missing(j)) {
        return(x[i, j, drop])
    } else if (!missing(i)) {
        return(x[i, , drop])
    } else {
        return(x[, j, drop])
    }
}

#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "sparseMatrix", i <- "nsparseVector", j <- "nsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "sparseMatrix", i <- "missing", j <- "nsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "sparseMatrix", i <- "nsparseVector", j <- "missing", drop="logical"), subset_generic_with_vector)

#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "sparseMatrix", i <- "nsparseVector", j <- "nsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "sparseMatrix", i <- "missing", j <- "nsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "sparseMatrix", i <- "nsparseVector", j <- "missing", drop="missing"), subset_generic_with_vector)



#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "matrix", i <- "nsparseVector", j <- "nsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "matrix", i <- "missing", j <- "nsparseVector", drop="logical"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "matrix", i <- "nsparseVector", j <- "missing", drop="logical"), subset_generic_with_vector)

#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "matrix", i <- "nsparseVector", j <- "nsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "matrix", i <- "missing", j <- "nsparseVector", drop="missing"), subset_generic_with_vector)
#' @rdname slice
#' @export
setMethod(`[`, signature(x <- "matrix", i <- "nsparseVector", j <- "missing", drop="missing"), subset_generic_with_vector)
