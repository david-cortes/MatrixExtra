library("testthat")
library("Matrix")
library("MatrixExtra")
context("Utility functions")

test_that("Removing zeros", {
    set.seed(1)
    X <- rsparsematrix(10, 5, .75, repr="T")
    X@x[5:8] <- 0
    Xt <- remove_sparse_zeros(X)
    Xr <- remove_sparse_zeros(as.csr.matrix(X))
    
    expect_equal(as.matrix(X), as.matrix(Xt))
    expect_equal(sum(Xt@x == 0), 0)
    
    expect_equal(as.matrix(X), unname(as.matrix(Xr)))
    expect_equal(sum(Xr@x == 0), 0)
    
    X@x[8:10] <- NA_real_
    Xt <- remove_sparse_zeros(X, na.rm=TRUE)
    Xr <- remove_sparse_zeros(as.csr.matrix(X), na.rm=TRUE)
    Xi <- unname(as.matrix(X))
    Xi[is.na(as.matrix(X))] <- 0
    
    expect_equal(Xi, as.matrix(Xt))
    expect_false(anyNA(Xt@x))
    
    expect_equal(Xi, unname(as.matrix(Xr)))
    expect_false(anyNA(Xr@x))
})

test_that("Sorting indices", {
    X <- new("dgRMatrix")
    X@p <- as.integer(c(0, 1, 4, 5, 6))
    X@j <- as.integer(c(4,  2,1,4,  1,  0))
    X@x <- c(-0.91, 0.14, -0.12, -0.12, 1.1, 0.66)
    X@Dim <- c(4L, 5L)
    
    X_copy <- deepcopy_sparse_object(X)
    indices <- X@j

    X_new <- sort_sparse_indices(X, copy=TRUE)
    expect_equal(X_new@j, as.integer(c(4,  1,2,4,  1,  0)))
    expect_equal(indices, as.integer(c(4,  2,1,4,  1,  0)))
    
    sort_sparse_indices(X, copy=FALSE)
    expect_equal(X@j, as.integer(c(4,  1,2,4,  1,  0)))
    expect_equal(X@j, indices)
})

test_that("Checking indices", {
    X <- new("dgRMatrix")
    X@p <- as.integer(c(0, 1, 4, 5, 6))
    X@j <- as.integer(c(4,  2,1,4,  1,  0))
    X@x <- c(-0.91, 0.14, -0.12, -0.12, 1.1, 0.66)
    X@Dim <- c(4L, 5L)
    
    X_copy <- X
    X <- check_sparse_matrix(X)
    expect_equal(X@j, as.integer(c(4,  1,2,4,  1,  0)))
    expect_equal(X_copy@j, as.integer(c(4,  2,1,4,  1,  0)))
    
    X@p <- as.integer(c(0, 1, 4, 5, 100))
    expect_error(check_sparse_matrix(X))
    X@p <- as.integer(c(0, 5, 4, 5, 6))
    expect_error(check_sparse_matrix(X))
    X@p <- as.integer(c(0, 1, NA, 5, 6))
    expect_error(check_sparse_matrix(X))
    X@p <- as.integer(c(0, -1, 4, 5, 6))
    expect_error(check_sparse_matrix(X))
    
    X@p <- as.integer(c(0, 1, 4, 5, 6))
    X@j <- as.integer(c(4,  1,2,4,  1,  10))
    expect_error(check_sparse_matrix(X))
    X@j <- as.integer(c(4,  1,2,4,  -1,  0))
    expect_error(check_sparse_matrix(X))
    X@j <- as.integer(c(4,  2,1,4,  1,  0))
    check_sparse_matrix(X)
})

test_that("Empty matrices", {
    X <- emptySparse(0, 1, format="R")
    expect_s4_class(X, "dgRMatrix")
    X <- emptySparse(1, 0, format="C", dtype="l")
    expect_s4_class(X, "lgCMatrix")
    X <- emptySparse(0, 0, format="T", dtype="n")
    expect_s4_class(X, "ngTMatrix")
    expect_error(suppressWarnings({X <- emptySparse(2^54, 1)}))
    expect_error({X <- emptySparse(format="Q")})
    expect_error({X <- emptySparse(dtype="i")})
})

test_that("Filter matrices", {
    set.seed(1)
    X <- rsparsematrix(nrow=20, ncol=10, density=0.3)
    Xcsr <- as.csr.matrix(X)
    Xcsc <- as.csc.matrix(X)
    Xcoo <- as.coo.matrix(X)
    svec <- as(Xcsr[1, ,drop=FALSE], "sparseVector")
    svec_num <- svec
    svec_num@i <- as.numeric(svec_num@i)
    
    X <- as.matrix(X)
    X[!((X == 0) | (X > 0.1))] <- 0
    
    res <- filterSparse(Xcsr, function(x) x > 0.1)
    expect_s4_class(res, "dgRMatrix")
    expect_equal(unname(X), unname(as.matrix(res)))
    
    res <- filterSparse(Xcsc, function(x) x > 0.1)
    expect_s4_class(res, "dgCMatrix")
    expect_equal(unname(X), unname(as.matrix(res)))
    
    res <- filterSparse(Xcoo, function(x) x > 0.1)
    expect_s4_class(res, "dgTMatrix")
    expect_equal(unname(X), unname(as.matrix(res)))
    
    res <- filterSparse(svec, function(x) x > 0.1)
    expect_s4_class(res, "dsparseVector")
    
    res <- filterSparse(svec_num, function(x) x > 0.1)
    expect_s4_class(res, "dsparseVector")
})
