#' @name MatrixExtra-options
#' @title MatrixExtra internal options
#' @description Controls some of the behaviors when calling method from
#' `MatrixExtra` which differ very significantly from `Matrix` and which
#' have the potential for breaking existing code in other packages.
#' 
#' The `Matrix` package has some particular choice of behaviors which can
#' inadvertently make operations very inefficient. Particularly:\itemize{
#' \item When transposing a sparse matrix in CSR or CSC formats, `Matrix` will
#' output a transpose in the same storage order. This is a slow operation, as it
#' requires duplicating the data, creating a new index, and sorting it. In general,
#' one might not care about the storage order of a sparse matrix until it comes the
#' time of performing an operation which is more efficient in one format, or one
#' might want to pass a sparse matrix object to a function from another package which
#' is not under one's control, for which a slow `t(.)` is undesirable.
#' 
#' `MatrixExtra` will instead make the default for `t(CSR)` and `t(CSC)` to output
#' in the opposite format (`t(CSR)` -> `CSC`; `t(CSC` -> `CSR`), which does not
#' involve any changes or duplication in the data, and it's thus much faster.
#' This in particular could break code in existing packages, particularly in
#' \bold{the `Matrix` package itself}.
#' The old `Matrix` behavior can be brought back through
#' `options("MatrixExtra.fast_transpose" = FALSE)`.
#' 
#' \item When selecting slices of a sparse matrix with only one row or only one
#' column, it will by default simplify them to a \bold{dense} vector, which goes
#' against the purpose of using sparse structures. With this being the default,
#' it's possible to inadvertently do this and end up losing all the benefits of
#' a sparse data structure.
#' 
#' `MatrixExtra` instead will make the default to simplify them to a \bold{sparse}
#' vector (see \link{slice}). It can be changed back to the `Matrix` behavior
#' through `options("MatrixExtra.drop_sparse" = FALSE)`.
#' 
#' \item When doing elementwise multiplication (`*`) of a sparse matrix by a
#' dense matrix or vice-versa, if the dense matrix has a `NA` or `NaN` element
#' in a coordinate at which the sparse matrix has no entry, `Matrix` will preserve
#' the `NA` in the resulting output, as does base R.
#' 
#' `MatrixExtra` will instead ignore such `NA`s, which makes the operation much
#' faster. The old `Matrix` behavior can be brought back through
#' `options("MatrixExtra.ignore_na" = FALSE)`.
#' 
#' \item When doing sparse-dense matrix multiplications, `Matrix` will run
#' single-threaded, which is typically undesirable as modern CPUs have the capacity
#' to perform such operations much faster by exploiting both multi-threading and SIMD
#' instructions.
#' 
#' `MatrixExtra` will by default use all the available threads in the system.
#' The number of threads can be controlled through `options("MatrixExtra.nthreads" = 1L)`.
#' }
#' 
#' These behaviors will change at the moment one loads `MatrixExtra` through
#' `library(MatrixExtra)`, save for the number of threads which will also be set
#' to the maximum when calling functions as e.g. `MatrixExtra::crossprod(x,y)`
#' 
#' The package provides two functions to make all these changes simultaneously:\itemize{
#' \item `restore_old_matrix_behavior`: Will change the transpose behavior to deep
#' transposes in the same format, drop behavior to dense vectors, preserve NAs that are
#' not in a sparse matrix in elementwise multiplication, and number of threads
#' to 1.
#' 
#' These all match with how the `Matrix` package behaves in such situations.
#' \item `set_new_matrix_behavior`: Will change the transpose behavior to shallow
#' transposes in the opposite format, drop behavior to sparse vectors, ignore NAs that
#' are not in a sparse matrix in elementwise multiplication, and number of
#' threads to the available number of threads in the system.
#' 
#' These all differ with how the `Matrix` package behaves in such situations.
#' }
NULL

#' @importFrom parallel detectCores

#' @rdname MatrixExtra-options
#' @export
set_new_matrix_behavior <- function() {
    options("MatrixExtra.drop_sparse" = TRUE)
    options("MatrixExtra.fast_transpose" = TRUE)
    options("MatrixExtra.ignore_na" = TRUE)
    options("MatrixExtra.nthreads" = parallel::detectCores())
}

#' @rdname MatrixExtra-options
#' @export
restore_old_matrix_behavior <- function() {
    options("MatrixExtra.drop_sparse" = FALSE)
    options("MatrixExtra.fast_transpose" = FALSE)
    options("MatrixExtra.ignore_na" = FALSE)
    options("MatrixExtra.nthreads" = 1L)
}


.onAttach <- function(libname, pkgname) {
    set_new_matrix_behavior()
    packageStartupMessage(
        "'MatrixExtra' modifies important behaviors from 'Matrix'. See ?MatrixExtra-options."
    )
}

.onLoad <- function(libname, pkgname) {
    restore_old_matrix_behavior()
    if (getOption("MatrixExtra.nthreads") == 1L)
        options("MatrixExtra.nthreads" = parallel::detectCores())
}
