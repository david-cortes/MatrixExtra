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
