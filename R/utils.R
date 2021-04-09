#' @title Sort the indices of a sparse matrix or sparse vector
#' @description Will sort the indices of a sparse matrix or sparse vector.
#'
#' In general, the indices of sparse CSR and CSC matrices are always meant to be sorted, and
#' it should be rare to have a matrix with unsorted indices when the matrix is the
#' output of some built-in operation from either `Matrix` or `MatrixExtra`, but when
#' the matrices are constructed manually this function can come in handy.
#'
#' \bold{Important:} the input matrix will be modified in-place, unless passing
#' `copy=TRUE`.
#' @param X A sparse matrix in CSR, CSC, or COO format; or a sparse vector
#' (from the `Matrix` package.)
#' @param copy Whether to make a deep copy of the indices and the values before sorting
#' them, so that the in-place modifications will not affect any potential external
#' references to the same arrays.
#' @param byrow When passing a COO matrix ("TsparseMatrix"), whether to sort it
#' with rows as the major axis, which differs from `Matrix` that usually
#' sorts them with columns as the major axis.
#' @return The same input `X` with its indices sorted, as invisible (no auto print).
#' Note that the input is itself modified, so there is no need to reassign it.
#' @export
sort_sparse_indices <- function(X, copy=FALSE, byrow=TRUE) {
    if (inherits(X, "RsparseMatrix")) {

        check_valid_matrix(X)
        if (copy) X@j <- deepcopy_int(X@j)

        if (inherits(X, "dsparseMatrix")) {
            if (copy) X@x <- deepcopy_num(X@x)
            sort_sparse_indices_numeric_known_ncol(
                X@p,
                X@j,
                X@x,
                ncol(X)
            )
        } else if (inherits(X, "lsparseMatrix")) {
            if (copy) X@x <- deepcopy_log(X@x)
            sort_sparse_indices_logical_known_ncol(
                X@p,
                X@j,
                X@x,
                ncol(X)
            )
        } else if (inherits(X, "nsparseMatrix")) {
            sort_sparse_indices_binary(
                X@p,
                X@j
            )
        } else {
            X <- as.csr.matrix(X)
            return(sort_sparse_indices(X, copy=copy))
        }

    } else if (inherits(X, "CsparseMatrix")) {

        check_valid_matrix(X)
        if (copy) X@i <- deepcopy_int(X@i)

        if (inherits(X, "dsparseMatrix")) {
            if (copy) X@x <- deepcopy_num(X@x)
            sort_sparse_indices_numeric_known_ncol(
                X@p,
                X@i,
                X@x,
                nrow(X)
            )
        } else if (inherits(X, "lsparseMatrix")) {
            if (copy) X@x <- deepcopy_log(X@x)
            sort_sparse_indices_logical_known_ncol(
                X@p,
                X@i,
                X@x,
                nrow(X)
            )
        } else if (inherits(X, "nsparseMatrix")) {
            sort_sparse_indices_binary(
                X@p,
                X@i
            )
        } else {
            X <- as.csc.matrix(X)
            return(sort_sparse_indices(X, copy=copy))
        }

    } else if (inherits(X, "TsparseMatrix")) {

        check_valid_matrix(X)
        
        if (!byrow)
            X <- t_shallow(X)
        
        if (copy) {
            X@i <- deepcopy_int(X@i)
            X@j <- deepcopy_int(X@j)
        }

        if (inherits(X, "dsparseMatrix")) {
            if (copy) X@x <- deepcopy_num(X@x)
            sort_coo_indices_numeric(
                X@i,
                X@j,
                X@x
            )
        } else if (inherits(X, "lsparseMatrix")) {
            if (copy) X@x <- deepcopy_log(X@x)
            sort_coo_indices_logical(
                X@i,
                X@j,
                X@x
            )
        } else if (inherits(X, "nsparseMatrix")) {
            sort_coo_indices_binary(
                X@i,
                X@j
            )
        } else {
            if (!byrow)
                X <- t_shallow(X)
            X <- as.coo.matrix(X)
            return(sort_sparse_indices(X, copy=copy, byrow=byrow))
        }
        
        if (!byrow)
            X <- t_shallow(X)

    } else if (inherits(X, "sparseVector")) {
        
        if (copy) X@i <- deepcopy_int(X@i)

        if (inherits(X, "dsparseVector")) {
            if (copy) X@x <- deepcopy_num(X@x)
            sort_vector_indices_numeric(
                X@i,
                X@x
            )
        } else if (inherits(X, "isparseVector")) {
            if (copy) X@x <- deepcopy_int(X@x)
            sort_vector_indices_integer(
                X@i,
                X@x
            )
        } else if (inherits(X, "lsparseVector")) {
            if (copy) X@x <- deepcopy_log(X@x)
            sort_vector_indices_logical(
                X@i,
                X@x
            )
        } else if (inherits(X, "nsparseVector")) {
            sort_vector_indices_binary(
                X@i
            )
        } else {
            X <- as(X, "dsparseVector")
            return(sort_sparse_indices(X, copy=copy))
        }

    } else {
        stop("Method is only applicable to sparse matrices in CSR, CSC, and COO formats, and to sparse vectors.")
    }
    return(invisible(X))
}

### TODO: complete this and add it where necessary
deepcopy_before_sort <- function(X, logical=FALSE, binary=FALSE) {
    if (logical)
        target <- c("lsparseMatrix", "lsparseVector")
    else
        target <- c("dsparseMatrix", "dsparseVector")
    
    if ((!inherits(X, target) ||
         (binary && !inherits(X, c("nsparseMatrix", "nsparseVector")))
         ) ||
        inherits(X, "symmetricMatrix") ||
        (.hasSlot(X, "diag") && X@diag != "N")
    ) {
        if (inherits(X, "RsparseMatrix")) {
            X@j <- deepcopy_int(X@j)
        } else if (inherits(X, "CsparseMatrix")) {
            X@i <- deepcopy_int(X@i)
        } else if (inherits(X, "TsparseMatrix")) {
            X@i <- deepcopy_int(X@i)
        }

        if (inherits(X, "dsparseMatrix")) {
            X@x <- deepcopy_num(X@x)
        } else if (inherits(X, "lsparseMatrix")) {
            X@x <- deepcopy_log(X@x)
        }
    }
    return(X)
}


#' @title Deep copy sparse matrices and vectors
#' @description Generates a deep copy of a sparse matrix or sparse vector
#' object, which can come useful when the matrix is to later be passed to
#' functions that will potentially modify it in-place, such as \link{sort_sparse_indices}.
#' @param X A sparse matrix or sparse vector from the `Matrix` package.
#' @return The same input `X` with the fields replaced with deep copies.
#' @export
deepcopy_sparse_object <- function(X) {
    if (inherits(X, "sparseVector")) {

        X@i <- deepcopy_int(X@i)
        X@length <- deepcopy_int(X@length)
        if (inherits(X, "dsparseVector")) {
            X@x <- deepcopy_num(X@x)
        } else if (inherits(X, "isparseVector")) {
            X@x <- deepcopy_int(X@x)
        } else if (inherits(X, "lsparseVector")) {
            X@x <- deepcopy_log(X@x)
        }

    } else if (inherits(X, "sparseMatrix")) {

        if (inherits(X, "TsparseMatrix")) {
            X@i <- deepcopy_int(X@p)
            X@j <- deepcopy_int(X@j)
        } else if (inherits(X, "CsparseMatrix")) {
            X@p <- deepcopy_int(X@p)
            X@i <- deepcopy_int(X@i)
        } else if (inherits(X, "RsparseMatrix")) {
            X@p <- deepcopy_int(X@p)
            X@j <- deepcopy_int(X@j)
        } else {
            throw_internal_error()
        }

        X@Dim <- deepcopy_int(X@Dim)
        if (inherits(X, "dsparseMatrix")) {
            X@x <- deepcopy_num(X@x)
        } else if (inherits(X, "lsparseMatrix")) {
            X@x <- deepcopy_log(X@x)
        }

        if (.hasSlot(X, "diag"))
            X@diag <- deepcopy_str(X@diag)
        if (.hasSlot(X, "uplo"))
            X@uplo <- deepcopy_str(X@uplo)

    } else {
        stop("Method is only applicable to sparse matrices and vectors.")
    }

    return(X)
}


#' @title Remove Zeros from a Sparse Matrix or Sparse Vector
#' @description Removes the entries in a sparse matrix or sparse vector which
#' have a value of zero but nevertheless are still among the object's values,
#' in any case there are any. Can also remove missing values if desired.
#' @param X A sparse matrix (COO, CSR, CSC) or sparse vector (any type)
#' from the `Matrix` package, whose values will be removed (left as non-present in the
#' sparse representation) if they are zeros.
#' @param na.rm Whether to also remove missing values (`NA` / `NaN`) from `X`.
#' @return The same matrix / vector X with its zeros removed from the sparse representation.
#' @export
remove_sparse_zeros <- function(X, na.rm=FALSE) {
    if (inherits(X, "sparseMatrix")) {

        if (inherits(X, "nsparseMatrix"))
            return(X)

        if (inherits(X, "TsparseMatrix")) {

            if (inherits(X, "dsparseMatrix")) {
                res <- remove_zero_valued_coo_numeric(X@i, X@j, X@x, na.rm)
            } else if (inherits(X, "lsparseMatrix")) {
                res <- remove_zero_valued_coo_logical(X@i, X@j, X@x, na.rm)
            } else {
                throw_internal_error()
            }
            X_attr <- attributes(X)
            X_attr$i <- res$ii
            X_attr$j <- res$jj
            X_attr$x <- res$xx
            attributes(X) <- X_attr
            return(X)

        } else if (inherits(X, "CsparseMatrix")) {

            if (inherits(X, "dsparseMatrix")) {
                res <- remove_zero_valued_csr_numeric(X@p, X@i, X@x, na.rm)
            } else if (inherits(X, "lsparseMatrix")) {
                res <- remove_zero_valued_csr_logical(X@p, X@i, X@x, na.rm)
            } else {
                throw_internal_error()
            }
            X_attr <- attributes(X)
            X_attr$p <- res$indptr
            X_attr$i <- res$indices
            X_attr$x <- res$values
            attributes(X) <- X_attr
            return(X)

        } else if (inherits(X, "RsparseMatrix")) {

            if (inherits(X, "dsparseMatrix")) {
                res <- remove_zero_valued_csr_numeric(X@p, X@j, X@x, na.rm)
            } else if (inherits(X, "lsparseMatrix")) {
                res <- remove_zero_valued_csr_logical(X@p, X@j, X@x, na.rm)
            } else {
                throw_internal_error()
            }
            X_attr <- attributes(X)
            X_attr$p <- res$indptr
            X_attr$j <- res$indices
            X_attr$x <- res$values
            attributes(X) <- X_attr
            return(X)

        } else {
            throw_internal_error()
        }

    } else if (inherits(X, "sparseVector")) {

        if (inherits(X, "nsparseVector"))
            return(X)

        if (inherits(X, "dsparseVector")) {
            res <- remove_zero_valued_svec_numeric(X@i, X@x, na.rm)
        } else if (inherits(X, "isparseVector")) {
            res <- remove_zero_valued_svec_integer(X@i, X@x, na.rm)
        } else if (inherits(X, "lsparseVector")) {
            res <- remove_zero_valued_svec_logical(X@i, X@x, na.rm)
        } else {
            throw_internal_error()
        }
        X_attr <- attributes(X)
        X_attr$i <- res$ii
        X_attr$x <- res$xx
        attributes(X) <- X_attr
        return(X)

    } else {
        stop("Function is only applicable to sparse matrices and sparse vectors.")
    }
}


check_valid_matrix <- function(X) {
    
    nrows <- nrow(X)
    if (is.na(nrows) || nrows < 0L)
        stop("Matrix has invalid number of rows.")
    ncols <- ncol(X)
    if (is.na(ncols) || ncols < 0L)
        stop("Matrix has invalid number of columns.")
    
    if (NROW(X@Dimnames[[1L]])) {
        if (length(X@Dimnames[[1L]]) != nrow(X))
            stop("Row names of matrix do not match with number of rows.")
    }
    if (NROW(X@Dimnames[[2L]])) {
        if (length(X@Dimnames[[2L]]) != ncol(X))
            stop("Column names of matrix do not match with number of columns.")
    }
    
    if (inherits(X, c("triangularMatrix", "symmetricMatrix"))) {
        if (nrow(X) != ncol(X))
            stop("Matrix type implies square shape, but number of rows and columns differ.")
        if (.hasSlot(X, "diag") && (NROW(X@diag) != 1L || !(X@diag %in% c("N", "U"))))
            stop("Matrix has invalid diagonal type.")
        if (NROW(X@uplo) != 1L || !(X@uplo %in% c("L", "U")))
            stop("Matrix has invalid 'uplo' field.")
    }

    if (inherits(X, "TsparseMatrix")) {
        
        if (length(X@i) != length(X@j))
            stop("Matrix is invalid (row and column indices have different length).")
        if (.hasSlot(X, "x") && length(X@x) != length(X@i))
            stop("Matrix is invalid (values and indices have different number of entries).")

    } else if (inherits(X, "RsparseMatrix")) {
        
        if (is.na(X@p[length(X@p)]))
            stop("Matrix is invalid (missing last index pointer, might indicate integer overflow).")
        if (.hasSlot(X, "x") && length(X@j) != length(X@x))
            stop("Matrix is invalid (lengths of indices and values differ).")
        if (length(X@p)-1L != X@Dim[1L])
            stop("Matrix is invalid ('p' doesn't match with dimension).")
        if (X@p[1L] != 0L || X@p[X@Dim[1L]+1L] != length(X@j))
            stop("Matrix is invalid ('p' has bad start/end.)")

    } else if (inherits(X, "CsparseMatrix")) {
        
        if (is.na(X@p[length(X@p)]))
            stop("Matrix is invalid (missing last index pointer, might indicate integer overflow).")
        if (.hasSlot(X, "x") && length(X@i) != length(X@x))
            stop("Matrix is invalid (lengths of indices and values differ).")
        if (length(X@p)-1L != X@Dim[2L])
            stop("Matrix is invalid ('p' doesn't match with dimension).")
        if (X@p[1L] != 0L || X@p[X@Dim[2L]+1L] != length(X@i))
            stop("Matrix is invalid ('p' has bad start/end.)")

    } else {
        
        throw_internal_error()

    }
}


#' @title Check if the underlying data behind a sparse matrix constitutes a valid object
#' @details Makes checks on the data contained in the sparse matrix or sparse vector
#' object for whether the data constitutes a valid matrix - e.g. indices must not be
#' negative or larger than the dimensions, index pointer must match with the indices,
#' etc.
#' 
#' As a short-hand, can also sort the matrix and remove zeros by calling the respective
#' functions \link{sort_sparse_indices} and \link{remove_sparse_zeros}.
#' 
#' A sparse matrix or sparse vector should never come out with invalid data from
#' functions from `Matrix` and `MatrixExtra`, with one exception: some `MatrixExtra`
#' functions might modify indices in-place, which can cause problems if the same
#' array of indices / values is also used by another matrix or R object elsewhere.
#' 
#' Otherwise, this function is aimed at making checks on matrices that are manually
#' constructed.
#' @param X A sparse matrix or sparse vector whose underlying arrays are to be checked.
#' @param sort Whether to sort the indices of `X` along the way. For this, will make
#' deep copies of the indices and values so that there are no issues with external
#' references being updated along the way.
#' @param remove_zeros Whether to remove entries in `X` which have a value of zero
#' but are nevertheless still present in the sparse representation.
#' @returns The same matrix or vector `X` (perhaps sorted or with zeros removed depending
#' on the parameters). If `X` contains data that doesn't make for a valid sparse matrix
#' (e.g. different number of values and indices), it will throw an error.
#' @export
check_sparse_matrix <- function(X, sort=TRUE, remove_zeros=TRUE) {
    
    if (inherits(X, "sparseMatrix")) {
        
        check_valid_matrix(X)
        
        if (inherits(X, "TsparseMatrix")) {
            res <- check_valid_coo_matrix(X@i, X@j, nrow(X), ncol(X))
        } else if (inherits(X, "RsparseMatrix")) {
            res <- check_valid_csr_matrix(X@p, X@j, nrow(X), ncol(X))
        } else if (inherits(X, "CsparseMatrix")) {
            res <- check_valid_csr_matrix(X@p, X@i, nrow(X), ncol(X))
        } else {
            throw_internal_error()
        }
        
        
    } else if (inherits(X, "sparseVector")) {

        if (is.na(X@length))
            stop("Vector has invalid length.")
        if (X@length < 0)
            stop("Vector has negative length.")
        
        if (!inherits(X, "nsparseVector")) {
            if (length(X@i) != length(X@x))
                stop("Vector indices and values have different length.")
        }

        res <- check_valid_svec(X@i, X@length)
        
    } else {
        stop("Function is only applicable to sparse matrices and sparse vectors.")
    }

    if (length(res)) {
        stop(res$err)
    }
    
    if (inherits(X, c("RsparseMatrix")))
        nnz_before <- length(X@j)
    else
        nnz_before <- length(X@i)
    if (remove_zeros) X <- remove_sparse_zeros(X)
    if (inherits(X, c("RsparseMatrix")))
        nnz_after <- length(X@j)
    else
        nnz_after <- length(X@i)
    if (sort) X <- sort_sparse_indices(X, copy=nnz_before == nnz_after)
    return(X)
}

throw_internal_error <- function () {
    stop("Internal error. Please open an issue in GitHub describing what you were doing.")
}

### TODO: one potential way of handling the issue of sorting in-place and rendering
### the inputs unusable, is by setting an attribute on the indices indicating which
### vectors, if any, are linked to it, and pre-check that modifying some indices
### would not leave an orphan vector of values elsewhere. These would have to be
### added everywhere the matrices and vectors are being used. For them, can use
### objects of class externalptr in R, so as to compare the pointer addresses,
### which also have the advantage of getting reset after begin deserialized in
### a new session.
### 
### Should also make an extra check on the matrix and make sure it gets sorted
### before being converted to a different dtype

### TODO: make a note of everything that could modify the inputs in-place,
### then fill in this

get_vector_pointer <- function(vector) {
    stop("Not yet implemented.")
    return(NULL)
}

add_linked_vec_to_indices <- function(indices, vector) {
    stop("Not yet implemented.")
    return(indices)
}

remove_linked_vec_from_indices <- function(indices, vector) {
    stop("Not yet implemented.")
    return(indices)
}

reset_linked_to_indices <- function(indices) {
    stop("Not yet implemented.")
    attributes(indices)$linked <- NULL
    return(indices)
}

can_modify_indices <- function(indices, vector=NULL) {
    stop("Not yet implemented.")
    return(TRUE)
}
