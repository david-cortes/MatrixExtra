library("testthat")
library("Matrix")
library("MatrixExtra")
context("Elementwise Matrix*svec")

add_NAs <- function(v, pct) {
    v[sample(length(v), floor(.1 * length(v)), replace=FALSE)] <- NA_real_
    return(v)
}

add_inf <- function(v, pct) {
    v[sample(length(v), floor(.1 * length(v)), replace=FALSE)] <- Inf
    v[sample(length(v), floor(.1 * length(v)), replace=FALSE)] <- (-Inf)
    return(v)
}

as.imatrix <- function(X) {
    X <- as.matrix(X)
    suppressWarnings(mode(X) <- "integer")
    return(X)
}

test_that("Exact shape", {
    set.seed(1)
    X <- rsparsematrix(100, 35, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(100, 1, .2))
    
    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
    
    expect_equal(unname(as.matrix(as.matrix(X) * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_equal(unname(as.matrix(as.imatrix(X) * v)),
                 unname(as.imatrix(X) * as.numeric(v)))
})

test_that("Recycled", {
    set.seed(1)
    X <- rsparsematrix(100, 35, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(25, 1, .2))

    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
    
    expect_equal(unname(as.matrix(as.matrix(X) * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_equal(unname(as.matrix(as.imatrix(X) * v)),
                 unname(as.imatrix(X) * as.numeric(v)))
})

test_that("Exact shape with NA and Inf", {
    set.seed(1)
    X <- rsparsematrix(10, 5, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(10, 1, .2))
    
    X@x <- add_NAs(X@x, .1)
    X@x <- add_inf(X@x, .1)
    v@x <- add_NAs(v@x, .1)
    v@x <- add_inf(v@x, .1)
    
    restore_old_matrix_behavior()
    
    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
    
    expect_equal(unname(as.matrix(as.matrix(X) * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_equal(unname(as.matrix(as.imatrix(X) * v)),
                 unname(as.imatrix(X) * as.numeric(v)))
})

test_that("Recycled with NA and Inf", {
    set.seed(1)
    X <- rsparsematrix(10, 5, .2, repr="R")
    v <- as.sparse.vector(rsparsematrix(10, 1, .2))
    
    X@x <- add_NAs(X@x, .1)
    X@x <- add_inf(X@x, .1)
    v@x <- add_NAs(v@x, .1)
    v@x <- add_inf(v@x, .1)

    restore_old_matrix_behavior()
    
    expect_equal(unname(as.matrix(X * v)),
                 unname(as.matrix(X) * as.numeric(v)))
    expect_s4_class(X * v, "dgRMatrix")
    
    expect_equal(unname(as.matrix(as.imatrix(X) * v)),
                 unname(as.imatrix(X) * as.numeric(v)))
})

test_that("Atypical recycles", {
    set.seed(1)
    X <- rsparsematrix(100, 35, .2, repr="R")
    v_factor_larger <- as.sparse.vector(rsparsematrix(200, 1, .2))
    v_uneven_smaller <- as.sparse.vector(rsparsematrix(30, 1, .2))
    v_uneven_larger <- as.sparse.vector(rsparsematrix(111, 1, .2))
    v_uneven_larger2 <- as.sparse.vector(rsparsematrix(222, 1, .2))
    v_full <- as.sparse.vector(rsparsematrix(nrow(X)*ncol(X), 1, .2))
    
    suppressWarnings({
        expect_equal(unname(as.matrix(as.matrix(X) * v_factor_larger)),
                     unname(as.matrix(X) * as.numeric(v_factor_larger)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_uneven_smaller)),
                     unname(as.matrix(X) * as.numeric(v_uneven_smaller)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_uneven_larger)),
                     unname(as.matrix(X) * as.numeric(v_uneven_larger)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_uneven_larger2)),
                     unname(as.matrix(X) * as.numeric(v_uneven_larger2)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_full)),
                     unname(as.matrix(X) * as.numeric(v_full)))
    })
    
    X@x <- add_NAs(X@x, .1)
    X@x <- add_inf(X@x, .1)
    v_factor_larger@x <- add_NAs(v_factor_larger@x, .1)
    v_factor_larger@x <- add_inf(v_factor_larger@x, .1)
    v_uneven_smaller@x <- add_NAs(v_uneven_smaller@x, .1)
    v_uneven_smaller@x <- add_inf(v_uneven_smaller@x, .1)
    v_uneven_larger@x <- add_NAs(v_uneven_larger@x, .1)
    v_uneven_larger@x <- add_inf(v_uneven_larger@x, .1)
    v_uneven_larger2@x <- add_NAs(v_uneven_larger2@x, .1)
    v_uneven_larger2@x <- add_inf(v_uneven_larger2@x, .1)
    v_full@x <- add_NAs(v_full@x, .1)
    v_full@x <- add_inf(v_full@x, .1)
    
    restore_old_matrix_behavior()
    suppressWarnings({
        expect_equal(unname(as.matrix(as.matrix(X) * v_factor_larger)),
                     unname(as.matrix(X) * as.numeric(v_factor_larger)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_uneven_smaller)),
                     unname(as.matrix(X) * as.numeric(v_uneven_smaller)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_uneven_larger)),
                     unname(as.matrix(X) * as.numeric(v_uneven_larger)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_uneven_larger2)),
                     unname(as.matrix(X) * as.numeric(v_uneven_larger2)))
        expect_equal(unname(as.matrix(as.matrix(X) * v_full)),
                     unname(as.matrix(X) * as.numeric(v_full)))
    })
})
