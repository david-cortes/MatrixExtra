# MatrixExtra

`MatrixExtra` is an R package which extends the sparse matrix and sparse vector types in the [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html) package, particularly the [CSR](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)) or `RsparseMatrix` formats (row-major), by providing methods that work natively and efficiently on them without converting them to another format along the way, such as slicing (selecting rows/columns) or concatenating by rows/columns, along with replacing some `Matrix` methods with more efficient versions, such as multi-threaded `<sparse, dense>` matrix multiplications, much faster slicing for all the sparse types, and faster elementwise addition/subtraction/multiplication, among others.

This package is based on code originally written for the [rsparse](https://github.com/rexyai/rsparse) package by Dmitriy Selivanov, and includes the MIT-licensed [robin-map](https://github.com/Tessil/robin-map) library within its source code.

# What's missing from Matrix

The `Matrix` package provides a rich set of sparse matrix and sparse vector classes with many methods and operators so that they could be used as drop-in replacements of base R's matrices. Unfortunately, the whole package is centered around the CSC format (`CsparseMatrix`, column-major), and calling methods and operators which in principle should be efficient in CSR or COO formats will imply first converting the whole matrix to CSC format (a slow and inefficient operation which duplicates the data), on which the operation might be less efficient due to the storage order.

*(Longer introduction in the [vignette](http://htmlpreview.github.io/?https://github.com/david-cortes/MatrixExtra/blob/master/inst/doc/Introducing_MatrixExtra.html))*

Examples:

```r
library(Matrix)
X = matrix(c(0,0,2,1, 0,3,0,0, 0,0,0,0), nrow=3, ncol=4, byrow=TRUE)
X = as(X, "RsparseMatrix")
X
```
```
3 x 4 sparse Matrix of class "dgRMatrix"
            
[1,] . . 2 1
[2,] . 3 . .
[3,] . . . .
```
```r
### This will forcibly convert the matrix to triplets
X[1:2, ]
```
```
2 x 4 sparse Matrix of class "dgTMatrix"
            
[1,] . . 2 1
[2,] . 3 . .
```
```r
### This will forcibly convert the matrix to CSC
rbind(X, X)
```
```
6 x 4 sparse Matrix of class "dgCMatrix"
            
[1,] . . 2 1
[2,] . 3 . .
[3,] . . . .
[4,] . . 2 1
[5,] . 3 . .
[6,] . . . .
```
```r
### This will forcibly convert the matrix to CSC
X * X
```
```
3 x 4 sparse Matrix of class "dgCMatrix"
            
[1,] . . 4 1
[2,] . 9 . .
[3,] . . . .
```

# Why is CSR needed?

The CSR sparse format is particularly useful when dealing with machine learning applications - e.g. splitting between a train and test set, tokenizing text features, multiplying a matrix by a vector of coefficients, calculating a gradient observation-by-observation, among others. Many stochastic optimization techniques and libraries (e.g. LibSVM, VowpalWabbit) require the inputs to be in CSR format or alike (see also [readsparse](https://www.github.com/david-cortes/readsparse)), which does not play well with the column-centric methods of `Matrix`.

In principle, one could stick with just the CSC format from `Matrix` and keep a mental map of the matrix as being transposed. This however gets complicated rather soon and is very prone to errors. Additionally, one might want to pass sparse matrices to another package whose code is outside of one's control, for which the storage format can make a large difference in performance.

# Installation

**Note:** This package greatly benefits from extra optimizations that aren't enabled by default for R packages. See [this guide](https://github.com/david-cortes/installing-optimized-libraries) for instructions on how to enable them.

```r
install.packages("MatrixExtra")
```

# Documentation

Documentation is available on [CRAN](https://cran.r-project.org/web/packages/MatrixExtra/index.html).

A package vignette is available [here](http://htmlpreview.github.io/?https://github.com/david-cortes/MatrixExtra/blob/master/inst/doc/Introducing_MatrixExtra.html) (and in CRAN).

# Features

* Multi-threaded matrix multiplications (`%*%`), `crossprod` and `tcrossprod` for many `<sparse,dense>` and `<dense,sparse>` types, including those from the [float](https://github.com/wrathematics/float) package.
* Slicing or subsetting (operator `[`) CSR matrices, along with faster slicing of CSC and COO, and slicing with sparse vectors.
* Efficient rbinding (concatenating by rows) and cbinding (concatenating by columns) for different sparse matrices and sparse vector types (e.g. `rbind(CSR, CSR)` and `cbind(CSR, CSR)`).
* Overloaded operators for `<RsparseMatrix, RsparseMatrix>`  and some `<RsparseMatrix, TsparseMatrix>` and `<sparse, dense>` types, such as `+`, `-`, `*`, `&`, `|`.
* Overloaded mathematical functions and operators which act only on the non-zero entries for CSR and COO matrices, such as `sqrt(CSR)` or `CSR * scalar`.
* Overloaded operators for `CSR`/`COO` and vectors, such as `CSR * vector`, `CSR ^ vector`, etc., which can mimick all the quirks of base R if needed.
* Convenience conversion functions between different sparse formats, and registered coercion methods between pairs which are not in `Matrix` (e.g. `matrix` -> `ngRMatrix` or `dgRMatrix` -> `lgCMatrix`).
* Fast transposes which work by outputting in the opposite format (CSR -> CSC and CSC -> CSR).
* Utility functions for sorting sparse indices and removing zero-valued entries.
* Convenience functions such as `mapSparse`, `emptySparse`, and a shorter `show` method for sparse objects.

# TODO

* Cover more cases in assignment operator (`[<-`).
* Matrix multiplications between `float32` and `sparseVector`.
* Better handling of dimension names of the output matrices.
* Outer products with sparse vectors.
* Try to port some parts to `Matrix`.
* Add timings against the methods from `Matrix`.
* Perhaps add methods specifically for symmetric and triangular matrices.

# Examples

```r
library(Matrix)
library(MatrixExtra)
X = matrix(c(0,0,2,1, 0,3,0,0, 0,0,0,0), nrow=3, ncol=4, byrow=TRUE)
X = as(X, "RsparseMatrix")
show(X)
options("MatrixExtra.quick_show" = FALSE)
show(X)

X[1:2, ]
X + X
X * X
X %*% 1:4
X * 1:3
rbind(X, X)
cbind(X, X)
sqrt(X)
diag(X)

as(as.matrix(X), "dgRMatrix")
as.csc.matrix(X)

### New and optional
set_new_matrix_behavior()
inherits(X[1,,drop=TRUE], "sparseVector")
inherits(t(X), "CsparseMatrix")
restore_old_matrix_behavior()
```
