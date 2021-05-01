#' @title Assignment operator for CSR matrices
#' @description Assign values to a CSR matrix.
#' Note: this will only be a relatively fast operation when
#' assigning contiguous row sequences. Only some of the potential
#' assignment cases to a CSR matrix are replaced here - for example,
#' cases that involve uneven recycling of vectors will be left to
#' the `Matrix` package.
#' @param x A CSR matrix whose values are to be replaced.
#' @param i The indices of the rows to replace.
#' @param j The indices of the columns to replace.
#' @param ... Not used
#' @param value The values to replace with.
#' @return The same `x` input with the values `[i,j]` set to `value`.
#' If the result is a full matrix (e.g. `x[,] <- 1`), the object will
#' be a dense matrix from base R.
#' @name assignment
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' set.seed(1)
#' X <- rsparsematrix(5, 3, .5, repr="R")
#' X[1:3] <- 0
#' X
NULL

throw_shape_err <- function() {
    stop("Values to assign does not match with matrix dimensions.")
}

assign_through_matrix <- function(x, i, j, value) {
    x <- as.coo.matrix(x)
    x[i, j] <- value
    x <- as.csr.matrix(x)
    return(x)
}

assign_csr_internal <- function(x, i, j, value, ij_properties=NULL) {

    if (!NROW(value) ||
        !inherits(value, c("numeric", "integer", "logical", "matrix",
                           "sparseMatrix", "sparseVector", "generalMatrix"))
    ) {
        stop("Invalid value to assign.")
    }
    if (inherits(value, "float32"))
        value <- float::dbl(value)
    if (inherits(value, "generalMatrix"))
        value <- as.matrix(value)
    if (is.matrix(value))
        value <- as.numeric(value)

    tot_x <- nrow(x) * ncol(x)


    if (is.null(ij_properties))
        ij_properties <- get_ij_properties(x, i ,j)

    i <- ij_properties$i
    j <- ij_properties$j
    all_i <- ij_properties$all_i
    all_j <- ij_properties$all_j
    i_is_seq <- ij_properties$i_is_seq
    j_is_seq <- ij_properties$j_is_seq
    i_is_rev_seq <- ij_properties$i_is_rev_seq
    j_is_rev_seq <- ij_properties$j_is_rev_seq
    n_row <- ij_properties$n_row
    n_col <- ij_properties$n_col

    if (!all_i || !i_is_seq || !i_is_rev_seq) {
        if (any(duplicated(i)))
            return(assign_through_matrix(x, i, j, value))
    }
    if (!all_j || !j_is_seq || !j_is_rev_seq) {
        if (any(duplicated(j)))
            return(assign_through_matrix(x, i, j, value))
    }

    if (inherits(value, c("numeric", "integer", "logical")) &&
        length(value) > 1 && length(i) > 1 && length(j) > 1 &&
        length(value) == length(i) * length(j)
    ) {
        value <- matrix(as.numeric(value), nrow=length(i), ncol=length(j))
    }
    if (inherits(value, "sparseVector") &&
        length(value) > 1 && length(i) > 1 && length(j) > 1 &&
        length(value) == (length(i) * length(j))
    ) {
        value <- Matrix(value, nrow=length(i), ncol=length(j), byrow=FALSE)
    }
    if (is.matrix(value) &&
        nrow(value) == length(i) && ncol(value) == length(j) &&
        length(i) > 1 && length(j) > 1 &&
        !all_i && !all_j
    ) {
        value <- as.csr.matrix(value)
    }


{
    if (is.vector(value) && length(value) == 1L)
    {
        if (!is.na(value) && value == 0)
        {
            if (all_i && all_j) {
                X_attr <- attributes(x)
                X_attr$x <- numeric()
                X_attr$j <- integer()
                X_attr$p <- integer(X_attr$Dim[1L]+1L)
                return(X_attr)
            }

            else if (all_j) {
                if (length(i) == 1L) {
                    res <- set_single_row_to_zero(x@p, x@j, x@x, i-1L)
                }

                else if (i_is_seq || i_is_rev_seq) {
                    imin <- ifelse(i_is_seq, i[1L], i[length(i)]) - 1L
                    imax <- ifelse(i_is_seq, i[length(i)], i[1L]) - 1L
                    res <- set_rowseq_to_zero(x@p, x@j, x@x, imin, imax)
                }

                else {
                    res <- set_arbitrary_rows_to_zero(x@p, x@j, x@x, i-1L)
                }
            }

            else if (all_i) {
                if (length(j) == 1L) {
                    res <- set_single_col_to_zero(x@p, x@j, x@x, j-1L)
                }

                else if (j_is_seq || j_is_rev_seq) {
                    jmin <- ifelse(j_is_seq, j[1L], j[length(j)]) - 1L
                    jmax <- ifelse(j_is_seq, j[length(j)], j[1L]) - 1L
                    res <- set_colseq_to_zero(x@p, x@j, x@x, jmin, jmax, ncol(x))
                }

                else {
                    res <- set_arbitrary_cols_to_zero(x@p, x@j, x@x, j-1L, ncol(x))
                }
            }

            else {
                if (length(i) == 1L && length(j) == 1L) {
                    res <- set_single_val_to_zero(x@p, x@j, x@x, i-1L, j-1L)
                }

                else if (length(j) == 1L) {
                    res <- set_arbitrary_rows_single_col_to_zero(x@p, x@j, x@x, i-1L, j-1L, ncol(x))
                }

                else if (length(i) == 1L) {
                    res <- set_single_row_arbitrary_cols_to_zero(x@p, x@j, x@x, i-1L, j-1L, ncol(x))
                }

                else {
                    res <- set_arbitrary_rows_arbitrary_cols_to_zero(x@p, x@j, x@x, i-1L, j-1L, ncol(x))
                }
            }
        }

        else ### value is not zero
        {
            value <- as.numeric(value)

            if (all_i && all_j) {
                warning("Warning: attempting to set all coordinates in a sparse matrix.")
                res <- matrix(value, nrow=nrow(x), ncol=ncol(x))
                if (!is.null(rownames(x)))
                    rownames(res) <- rownames(x)
                if (!is.null(colnames(x)))
                    colnames(res) <- colnames(x)
                return(res)
            }

            else if (all_j) {
                if (length(i) == 1L) {
                    res <- set_single_row_to_const(x@p, x@j, x@x, ncol(x), i-1L, value)
                }

                else if (i_is_seq || i_is_rev_seq) {
                    imin <- ifelse(i_is_seq, i[1L], i[length(i)]) - 1L
                    imax <- ifelse(i_is_seq, i[length(i)], i[1L]) - 1L
                    res <- set_rowseq_to_const(x@p, x@j, x@x, imin, imax, ncol(x), value)
                }

                else {
                    res <- set_arbitrary_rows_to_const(x@p, x@j, x@x, i-1L, ncol(x), value)
                }
            }

            else if (all_i) {
                if (length(j) == 1L) {
                    res <- set_single_col_to_const(x@p, x@j, x@x, ncol(x), j-1L, value)
                }

                else if (j_is_seq || j_is_rev_seq) {
                    jmin <- ifelse(j_is_seq, j[1L], j[length(j)]) - 1L
                    jmax <- ifelse(j_is_seq, j[length(j)], j[1L]) - 1L
                    res <- set_colseq_to_const(x@p, x@j, x@x, jmin, jmax, ncol(x), value)
                }

                else {
                    res <- set_arbitrary_cols_to_const(x@p, x@j, x@x, j-1L, ncol(x), value)
                }
            }

            else {
                if (length(i) == 1L && length(j) == 1L) {
                    res <- set_single_val_to_const(x@p, x@j, x@x, ncol(x), i-1L, j-1L, value)
                }

                else if (length(j) == 1L) {
                    res <- set_arbitrary_rows_single_col_to_const(x@p, x@j, x@x, i-1L, j-1L, value, ncol(x))
                }

                else if (length(i) == 1L) {
                    res <- set_single_row_arbitrary_cols_to_const(x@p, x@j, x@x, i-1L, j-1L, ncol(x), value)
                }

                else {
                    res <- set_arbitrary_rows_arbitrary_cols_to_const(x@p, x@j, x@x, i-1L, j-1L, ncol(x), value)
                }
            }
        }
    }

    else if (is.vector(value))
    {
        value <- as.numeric(value)

        if (all_i && all_j) {
            if ((length(value) > tot_x) || (tot_x %% length(value)) != 0)
                throw_shape_err()
            warning("Warning: attempting to set all coordinates of a sparse matrix.")
            res <- matrix(value, nrow=nrow(x), ncol=ncol(x))
            if (!is.null(rownames(x)))
                rownames(res) <- rownames(x)
            if (!is.null(colnames(x)))
                colnames(res) <- colnames(x)
            return(res)
        }

        else if (all_j) {

            if (length(i) == 1L) {
                if (length(value) > ncol(x) || (ncol(x) %% length(value)) != 0)
                    throw_shape_err()
                res <- set_single_row_to_rowvec(x@p, x@j, x@x, ncol(x), i-1L, value)
            }

            else {
                return(assign_through_matrix(x, i, j, value))
            }

        }

        else if (all_i) {

            if (length(j) == 1L) {
                if (length(value) > nrow(x) || (nrow(x) %% length(value)) != 0)
                    throw_shape_err()
                res <- set_single_col_to_colvec(x@p, x@j, x@x, ncol(x), j-1L, value)
            }

            else {
                return(assign_through_matrix(x, i, j, value))
            }

        }

        else {
            if (length(i) * length(j) != length(value))
                throw_shape_err()

            return(assign_through_matrix(x, i, j, value))
        }
    }

    else if (inherits(value, "sparseVector"))
    {
        if (length(value) == 0L)
            stop("Invalid value to assign to matrix.")
        if (length(value) == 1L)
            return(assign_csr_internal(x, i, j, as.numeric(value), ij_properties))

        if (all_i && all_j) {
            if ((length(value) > tot_x) || (tot_x %% length(value)) != 0)
                throw_shape_err()
            return(assign_through_matrix(x, i, j, value))
        }

        else if (all_j) {
            if (length(i) == 1L) {
                if (length(value) > ncol(x) || (ncol(x) %% length(value)) != 0)
                    throw_shape_err()
                if (length(value@i) == 0L)
                    return(assign_csr_internal(x, i, j, 0, ij_properties))
                res <- set_single_row_to_svec(x@p, x@j, x@x, ncol(x), i-1L, value@i-1L, value@x, length(value))
            }

            else {
                return(assign_through_matrix(x, i, j, value))
            }
        }

        else if (all_i) {
            if (length(j) == 1L) {
                if (length(value) > nrow(x) || (nrow(x) %% length(value)) != 0)
                    throw_shape_err()
                if (length(value@i) == 0L)
                    return(assign_csr_internal(x, i, j, 0, ij_properties))

                inplace_sort <- getOption("MatrixExtra.inplace_sort", default=FALSE)
                if (inplace_sort)
                    value <- deepcopy_before_sort(value)
                value <- sort_sparse_indices(value, copy=!inplace_sort)

                res <- set_single_col_to_svec(x@p, x@j, x@x, ncol(x), j-1L, value@i-1L, value@x, length(value))
            }

            else {
                return(assign_through_matrix(x, i, j, value))
            }
        }

        else {
            return(assign_through_matrix(x, i, j, value))
        }
    }

    else if (inherits(value, "sparseMatrix"))
    {
        tot_val <- nrow(value) * ncol(value)
        if (length(value@x) == 0L)
            return(assign_csr_internal(x, i, j, 0, ij_properties))
        
        if (all_i && all_j) {
            if (tot_val > tot_x || (tot_x %% tot_val) != 0)
                throw_shape_err()
            if (nrow(x) == nrow(value) && ncol(x) == ncol(value))
                return(as.csr.matrix(value))
            value <- as.sparse.vector(value)
            return(assign_csr_internal(x, i, j, value, ij_properties))
        }

        else if (all_j) {
            if (length(i) == 1L) {
                if (tot_val > ncol(x) || (ncol(x) %% tot_val) != 0)
                    throw_shape_err()
                if (nrow(value) != 1L && ncol(value) != 1L) {
                    value <- as.sparse.vector(value)
                    return(assign_csr_internal(x, i, j, value, ij_properties))
                }
                if (nrow(value) == 1L && ncol(value) == 1L) {
                    value <- as.numeric(as.matrix(value))
                    return(assign_csr_internal(x, i, j, value, ij_properties))
                }
                
                if (nrow(value) == 1L) {
                    if (inherits(value, "CsparseMatrix")) {
                        value <- as.sparse.vector(value)
                        return(assign_csr_internal(x, i, j, value, ij_properties))
                    }
                    if (inherits(value, c("RsparseMatrix", "TsparseMatrix"))) {
                        res <- set_single_row_to_svec(x@p, x@j, x@x, ncol(x), i-1L, value@j, value@x, ncol(value))
                    }
                    else {
                        throw_internal_error()
                    }
                }

                else if (ncol(value) == 1L) {
                    if (inherits(value, "RsparseMatrix")) {
                        value <- as.sparse.vector(value)
                        return(assign_csr_internal(x, i, j, value, ij_properties))
                    }
                    if (inherits(value, c("CsparseMatrix", "TsparseMatrix"))) {
                        res <- set_single_row_to_svec(x@p, x@j, x@x, ncol(x), i-1L, value@i, value@x, ncol(value))
                    }
                    else {
                        throw_internal_error()
                    }
                }

                else {
                    throw_internal_error()
                }
            }

            else if ((i_is_seq || i_is_rev_seq)) {

                if (nrow(value) == length(i) && ncol(value) == ncol(x)) {
                    value <- as.csr.matrix(value)
                    imin <- ifelse(i_is_seq, i[1L], i[length(i)])
                    imax <- ifelse(i_is_seq, i[length(i)], i[1L])
                    if (i_is_rev_seq)
                        value <- value[seq(length(i), 1L), , drop=FALSE]
                    res <- set_rowseq_to_smat(x@p, x@j, x@x, imin-1L, imax-1L, value@p, value@j, value@x)
                }

                else {
                    value <- as.sparse.vector(value)
                    return(assign_csr_internal(x, i, j, value, ij_properties))
                }
            }

            else {
                value <- as.csr.matrix(value)
                if (length(i) == nrow(x))
                    return(value[order(i), , drop=FALSE])
                if (!check_is_sorted(i)) {
                    argsorted <- order(i)
                    value <- value[argsorted, , drop=FALSE]
                    i <- i[argsorted]
                }

                res <- set_arbitrary_rows_to_smat(x@p, x@j, x@x, i-1L, value@p, value@j, value@x)
            }
        }

        else if (all_i) {
            if (length(j) == 1L) {
                if (tot_val > nrow(x) || (nrow(x) %% tot_val) != 0)
                    throw_shape_err()
                if (nrow(value) != 1L && ncol(value) != 1L) {
                    value <- as.sparse.vector(value)
                    return(assign_csr_internal(x, i, j, value, ij_properties))
                }
                if (nrow(value) == 1L && ncol(value) == 1L) {
                    value <- as.numeric(as.matrix(value))
                    return(assign_csr_internal(x, i, j, value, ij_properties))
                }

                if (nrow(value) == 1L) {
                    if (inherits(value, "CsparseMatrix")) {
                        value <- as.sparse.vector(value)
                        return(assign_csr_internal(x, i, j, value, ij_properties))
                    }
                    if (inherits(value, c("RsparseMatrix", "TsparseMatrix"))) {
                        res <- set_single_col_to_svec(x@p, x@j, x@x, ncol(x), j-1L, value@j, value@x, ncol(value))
                    }
                    else {
                        throw_internal_error()
                    }
                }

                else if (ncol(value) == 1L) {
                    if (inherits(value, "RsparseMatrix")) {
                        value <- as.sparse.vector(value)
                        return(assign_csr_internal(x, i, j, value, ij_properties))
                    }
                    if (inherits(value, c("CsparseMatrix", "TsparseMatrix"))) {
                        res <- set_single_col_to_svec(x@p, x@j, x@x, ncol(x), j-1L, value@i, value@x, ncol(value))
                    }
                    else {
                        throw_internal_error()
                    }
                }

                else {
                    throw_internal_error()
                }
            }

            else {
                return(assign_through_matrix(x, i, j, value))
            }
        }

        else {
            return(assign_through_matrix(x, i, j, value))
        }
    }

    else
    {
        if (!is.null(dim(value)))
            value <- as.matrix(value)
        else
            value <- as.numeric(value)
        return(assign_csr_internal(x, i, j, value, ij_properties))
    }

}

    X_attr <- attributes(x)
    X_attr$p <- res$indptr
    X_attr$j <- res$indices
    X_attr$x <- res$values
    attributes(x) <- X_attr
    return(x)
}

assign_csr <- function(x, i, j, ..., value) {
    if (length(list(...)))
        warning("Warning: unused arguments in assignment operator.")
    return(assign_csr_internal(x=x, i=i, j=j, value=value))
}

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="index", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="index", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="missing", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="missing", value="replValue"), assign_csr)



#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="index", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="index", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="missing", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="missing", value="sparseVector"), assign_csr)



#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="nsparseVector", j="index", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="nsparseVector", j="missing", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="nsparseVector", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="nsparseVector", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="nsparseVector", j="nsparseVector", value="replValue"), assign_csr)


#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="nsparseVector", j="index", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="nsparseVector", j="missing", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="nsparseVector", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="nsparseVector", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="nsparseVector", j="nsparseVector", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseVector", j="index", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseVector", j="missing", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="lsparseVector", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="lsparseVector", value="replValue"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseVector", j="lsparseVector", value="replValue"), assign_csr)


#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseVector", j="index", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseVector", j="missing", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="lsparseVector", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="lsparseVector", value="sparseVector"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseVector", j="lsparseVector", value="sparseVector"), assign_csr)


#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseMatrix", j="index", value="sparseMatrix"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseMatrix", j="missing", value="sparseMatrix"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="index", j="lsparseMatrix", value="sparseMatrix"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="missing", j="lsparseMatrix", value="sparseMatrix"), assign_csr)

#' @rdname assignment
#' @export
setMethod("[<-", signature(x="dgRMatrix", i="lsparseMatrix", j="lsparseMatrix", value="sparseMatrix"), assign_csr)
