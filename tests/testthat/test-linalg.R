library("testthat")
library("Matrix")
library("MatrixExtra")
context("linalg")

test_that("norm", {
    set.seed(1)
    X <- matrix(rnorm(100*5), nrow=5L)
    Xcsr <- as.csr.matrix(X)
    Xcsc <- as.csc.matrix(X)
    
    expect_equal(norm(Xcsr, "O"), norm(Xcsc, "O"))
    expect_equal(norm(Xcsr, "o"), norm(Xcsc, "o"))
    expect_equal(norm(Xcsr, "1"), norm(Xcsc, "1"))
    expect_equal(norm(Xcsr, "I"), norm(Xcsc, "I"))
    expect_equal(norm(Xcsr, "i"), norm(Xcsc, "i"))
    expect_equal(norm(Xcsr, "F"), norm(Xcsc, "F"))
    expect_equal(norm(Xcsr, "f"), norm(Xcsc, "f"))
    expect_equal(norm(Xcsr, "2"), norm(Xcsc, "2"))
})

test_that("diag", {
    set.seed(1)
    v <- rnorm(100)
    Xsquare <- matrix(v, nrow=10)
    Xwide <- matrix(v, nrow=5)
    Xdeep <- matrix(v, ncol=5)
    
    expect_equal(diag(as.csr.matrix(Xsquare)), diag(Xsquare))
    expect_equal(diag(as.csr.matrix(Xwide)), diag(Xwide))
    expect_equal(diag(as.csr.matrix(Xdeep)), diag(Xdeep))
})

test_that("diag assign", {
    set.seed(1)
    v <- rnorm(100)
    Xsquare <- matrix(v, nrow=10)
    Xwide <- matrix(v, nrow=5)
    Xdeep <- matrix(v, ncol=5)
    
    test_assign_diag <- function(X, v) {
        Xdense <- as.matrix(X)
        diag(Xdense) <- v
        diag(X) <- v
        expect_equal(unname(as.matrix(X)), unname(X))
    }
    
    test_assign_diag(Xsquare, 1)
    test_assign_diag(Xsquare, NA)
    test_assign_diag(Xsquare, 20*diag(Xsquare))
    expect_error(diag(Xsquare) <- 1:2)
    
    test_assign_diag(Xwide, 1)
    test_assign_diag(Xwide, NA)
    test_assign_diag(Xwide, 20*diag(Xwide))
    expect_error(diag(Xwide) <- 1:2)
    
    test_assign_diag(Xdeep, 1)
    test_assign_diag(Xdeep, NA)
    test_assign_diag(Xdeep, 20*diag(Xdeep))
    expect_error(diag(Xdeep) <- 1:2)
})
