#' @title Sort the indices of a sparse matrix or sparse vector
#' @description Will sort the indices of a sparse matrix or sparse vector.
#'
#' In general, the indices of sparse CSR and CSC matrices are always meant to be sorted, and
#' it should be rare to have a matrix with unsorted indices when the matrix is the
#' output of some built-in operation from either `Matrix` or `MatrixExtra`, but when
#' the matrices are constructed manually this function can come in handy.
#'
#' \bold{Important:} the input matrix will be modified in-place.
#' @details When passing a COO matrix ("TsparseMatrix"), it will sort it
#' with rows as the major axis, which differs from `Matrix` that usually
#' sorts them with columns as the major axis.
#' @param X A sparse matrix in CSR, CSC, or COO format; or a sparse vector
#' (from the `Matrix` package.)
#' @return The same input `X` with its indices sorted, as invisible (no auto print).
#' Note that the input is itself modified, so there is no need to reassign it.
#' @export
sort_sparse_indices <- function(X) {
    if (inherits(X, "RsparseMatrix")) {

        check_valid_matrix(X)

        if (inherits(X, "dsparseMatrix")) {
            sort_sparse_indices_numeric_known_ncol(
                X@p,
                X@j,
                X@x,
                ncol(X)
            )
        } else if (inherits(X, "lsparseMatrix")) {
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
            return(sort_sparse_indices(X))
        }

    } else if (inherits(X, "CsparseMatrix")) {

        check_valid_matrix(X)

        if (inherits(X, "dsparseMatrix")) {
            sort_sparse_indices_numeric_known_ncol(
                X@p,
                X@i,
                X@x,
                nrow(X)
            )
        } else if (inherits(X, "lsparseMatrix")) {
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
            return(sort_sparse_indices(X))
        }

    } else if (inherits(X, "TsparseMatrix")) {

        check_valid_matrix(X)

        if (inherits(X, "dsparseMatrix")) {
            sort_coo_indices_numeric(
                X@i,
                X@j,
                X@x
            )
        } else if (inherits(X, "lsparseMatrix")) {
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
            X <- as.coo.matrix(X)
            return(sort_sparse_indices(X))
        }

    } else if (inherits(X, "sparseVector")) {

        if (inherits(X, "dsparseVector")) {
            sort_vector_indices_numeric(
                X@i,
                X@x
            )
        } else if (inherits(X, "isparseVector")) {
            sort_vector_indices_integer(
                X@i,
                X@x
            )
        } else if (inherits(X, "lsparseVector")) {
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
            return(sort_sparse_indices(X))
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

### TODO: add tests for this
#' @title Deep copy sparse matrices and vectors
#' @description Generates a deep copy of a sparse matrix or sparse vector
#' object, which can come useful when the matrix is to later be passed to
#' functions that will potentially modify it in-place, such as \link{sort_sparse_indices}.
#' @param X A sparse matrix or sparse vector from the `Matrix` package.
#' @return The same input `X` with the fields replaced with deep copies.
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

### TODO: maybe make this an exportable function with an 'extensive check' parameter

check_valid_matrix <- function(X) {
    
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
        
        if (.hasSlot(X, "x") && length(X@j) != length(X@x))
            stop("Matrix is invalid (lengths of indices and values differ).")
        if (length(X@p)-1L != X@Dim[1L])
            stop("Matrix is invalid ('p' doesn't match with dimension).")
        if (X@p[1L] != 0L || X@p[X@Dim[1L]+1L] != length(X@j))
            stop("Matrix is invalid ('p' has bad start/end.)")

    } else if (inherits(X, "CsparseMatrix")) {
        
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

