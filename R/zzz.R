#' @name MatrixExtra-options
#' @title MatrixExtra internal options
#' @description Controls some of the behaviors when calling methods from
#' `MatrixExtra` which differ very significantly from the same methods in
#' `Matrix` and which have the potential for breaking existing code in other packages.
#' 
#' For example, the function `Matrix::sparse.model.matrix` might throw errors
#' if it is called after loading `library(MatrixExtra)` with the default options,
#' and the default transpose behavior needs to be modified from `MatrixExtra` to
#' make it work again.
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
#' For example, the function `Matrix::sparse.model.matrix` is unlikely to work
#' under the default fast transpose option, and getting it to work again will
#' require setting `options("MatrixExtra.fast_transpose" = FALSE)` or calling
#' `MatrixExtra::restore_old_matrix_behavior()`.
#' \item When selecting slices of a sparse matrix with only one row or only one
#' column, `Matrix` will by default simplify them to a \bold{dense} vector, which goes
#' against the purpose of using sparse structures. With this being the default,
#' it's possible to inadvertently do this and end up losing all the benefits of
#' a sparse data structure.
#' 
#' `MatrixExtra` instead will make the default to simplify them to a \bold{sparse}
#' vector (see \link{slice}). It can be changed back to the `Matrix` behavior
#' through `options("MatrixExtra.drop_sparse" = FALSE)`.
#' 
#' \item When doing some operations which require the indices of sparse matrices to
#' be sorted, such as elementwise addition and multiplication, `Matrix` will create
#' deep copies of the indices and values, which takes extra time but ensures that
#' no external references to the same objects are impacted, thus playing well with
#' R's pass-by-value and copy-on-modify semantics, but at the expense of extra time
#' and memory requirements.
#' 
#' `MatrixExtra` will instead sort the indices in-place when needed, only doing
#' deep copies when the values are to be cast to a different type, so that the
#' object being passed to a given function or operator will always remain usable
#' (having the same values at the same indices). However, when the indices of a given
#' object are unsorted and get sorted in-place, if any other R object has references
#' to the same vector of indices separately from the vector of values or vice-versa,
#' those objects will see a new vector which is modified, potentially rendering
#' those external objects unusable.
#' 
#' If this behavior is problematic, it can be changed back to always making deep copies
#' through `options("MatrixExtra.inplace_sort" = FALSE)`.
#' 
#' \item When doing elementwise multiplication (`*`) of a sparse matrix by a
#' dense matrix or vice-versa, if the dense matrix has a `NA` or `NaN` element
#' in a coordinate at which the sparse matrix has no entry, `Matrix` will preserve
#' the `NA` in the resulting output, as does base R. This also applies with the
#' propagation of 1s when using the `^` operator with zeros in the RHS, and with
#' propagation of infinities in division by zero.
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
#' 
#' \item When calling method `show` on a sparse matrix object (for example, by typing the
#' corresponding variable name in an R console and pressing 'Enter', `Matrix` will
#' take a glance at the object by printing the values in some or all of the rows and columns,
#' typically involving a long output which, for a mid-to-large sparse matrix, will typically be
#' larger than the available screen space and can make it very inconvenient to quickly inspect
#' sparse objects.
#' 
#' `MatrixExtra` will instead override this method with a quicker one that will only show
#' the type, dimensions, and number of entries in a sparse object, without printing the
#' values per-se as `Matrix` would do. Note that this only overrides the `show` method, while
#' `print` remains the same as it was. It can be changed back to the same `show` that is
#' provided by `Matrix` through `options("MatrixExtra.quick_show" = FALSE)`.
#' }
#' 
#' These behaviors will change at the moment one loads `MatrixExtra` through
#' `library(MatrixExtra)`, save for the number of threads which will also be set
#' to the maximum when calling functions as e.g. `MatrixExtra::<function>(<args>)`
#' 
#' The package provides two short-hand functions to switch all these options
#' simultaneously:\itemize{
#' \item `restore_old_matrix_behavior`: Will change the transpose behavior to deep
#' transposes in the same format, drop behavior to dense vectors, avoid in-place
#' sorting of sparse indices, preserve NAs that are not in a sparse matrix in
#' elementwise multiplication, and set the number of threads to 1.
#' 
#' These all match with how the `Matrix` package behaves in such situations.
#' \item `set_new_matrix_behavior`: Will change the transpose behavior to shallow
#' transposes in the opposite format, drop behavior to sparse vectors, sort sparse
#' indices in-place when possible, ignore NAs that are not in a sparse matrix in
#' elementwise multiplication, and set the number of threads to the available
#' number of threads in the system.
#' 
#' These all differ with how the `Matrix` package behaves in such situations.
#' }
#' @section All options:
#' List of options with the short-hand command to make them match with `Matrix`.
#' \itemize{
#' \item `options("MatrixExtra.fast_transpose" = FALSE)` :
#' Option for behavior of `t(CSR)` and `t(CSC)`.
#' \item `options("MatrixExtra.drop_sparse" = FALSE)` :
#' Option for behavior of `drop=TRUE` when subsetting.
#' \item `options("MatrixExtra.inplace_sort" = FALSE)` :
#' Option for behavior of inplace-vs-copy operations.
#' \item `options("MatrixExtra.ignore_na" = FALSE)` :
#' Option for behavior of `NA`s in elementwise
#' sparse-dense multipliction.
#' \item `options("MatrixExtra.nthreads" = 1L)` :
#' Number of parallel threads to use in sparse-dense matrix
#' \item `options("MatrixExtra.quick_show" = FALSE)` :
#' Option for behavior of `show` method for sparse objects.
#' }
#' @return No return value, called for side effects.
NULL

#' @importFrom parallel detectCores

#' @rdname MatrixExtra-options
#' @export
set_new_matrix_behavior <- function() {
    options("MatrixExtra.drop_sparse" = TRUE)
    options("MatrixExtra.fast_transpose" = TRUE)
    options("MatrixExtra.inplace_sort" = TRUE)
    options("MatrixExtra.ignore_na" = TRUE)
    options("MatrixExtra.nthreads" = parallel::detectCores())
    options("MatrixExtra.quick_show" = TRUE)
}

#' @rdname MatrixExtra-options
#' @export
restore_old_matrix_behavior <- function() {
    options("MatrixExtra.drop_sparse" = FALSE)
    options("MatrixExtra.fast_transpose" = FALSE)
    options("MatrixExtra.inplace_sort" = FALSE)
    options("MatrixExtra.ignore_na" = FALSE)
    options("MatrixExtra.nthreads" = 1L)
    options("MatrixExtra.quick_show" = FALSE)
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
