library("testthat")
library("Matrix")
library("MatrixExtra")
context("Elementwise CSR*svec")

test_that("CSR and Exact shape", {
    set.seed(1)
    X <- rsparsematrix(100, 35, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(100, 1, .2))
    
    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
})

test_that("CSR and Recycled", {
    set.seed(1)
    X <- rsparsematrix(100, 35, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(25, 1, .2))

    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
})

test_that("CSR and Exact shape with NA and Inf", {
    set.seed(1)
    X <- rsparsematrix(10, 5, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(10, 1, .2))
    
    nnz_X <- length(X@x)
    X@x[sample(nnz_X, floor(.1 * nnz_X), replace=FALSE)] <- NA_real_
    X@x[sample(nnz_X, floor(.1 * nnz_X), replace=FALSE)] <- Inf
    X@x[sample(nnz_X, floor(.1 * nnz_X), replace=FALSE)] <- (-Inf)
    
    nnz_v <- length(v@x)
    v@x[sample(nnz_v, floor(.1 * nnz_v), replace=FALSE)] <- NA_real_
    v@x[sample(nnz_v, floor(.1 * nnz_v), replace=FALSE)] <- Inf
    v@x[sample(nnz_v, floor(.1 * nnz_v), replace=FALSE)] <- (-Inf)
    
    restore_old_matrix_behavior()
    
    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
})

test_that("CSR and Recycled with NA and Inf", {
    set.seed(1)
    X <- rsparsematrix(10, 5, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(10, 1, .2))
    
    nnz_X <- length(X@x)
    X@x[sample(nnz_X, floor(.1 * nnz_X), replace=FALSE)] <- NA_real_
    X@x[sample(nnz_X, floor(.1 * nnz_X), replace=FALSE)] <- Inf
    X@x[sample(nnz_X, floor(.1 * nnz_X), replace=FALSE)] <- (-Inf)
    
    nnz_v <- length(v@x)
    v@x[sample(nnz_v, floor(.1 * nnz_v), replace=FALSE)] <- NA_real_
    v@x[sample(nnz_v, floor(.1 * nnz_v), replace=FALSE)] <- Inf
    v@x[sample(nnz_v, floor(.1 * nnz_v), replace=FALSE)] <- (-Inf)
    
    restore_old_matrix_behavior()
    
    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
})

