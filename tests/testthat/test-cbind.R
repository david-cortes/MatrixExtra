library("testthat")
library("Matrix")
library("MatrixExtra")
context("cbinding")

test_that("cbind COO", {
    set.seed(1)
    X <- as.coo.matrix(rsparsematrix(100, 50, .3))
    Y <- as.coo.matrix(rsparsematrix(100, 20, .3))

    expect_s4_class(cbind(X, Y), "dgTMatrix")
    expect_equal(as.matrix(cbind(X, Y)),
                 cbind(as.matrix(X), as.matrix(Y)))
    expect_equal(unname(as.matrix(cbind(as.csr.matrix(X), Y))),
                 cbind(as.matrix(X), as.matrix(Y)))
    expect_equal(unname(as.matrix(cbind(X, as.csr.matrix(Y)))),
                 cbind(as.matrix(X), as.matrix(Y)))
    expect_equal(unname(as.matrix(cbind(as.csc.matrix(X), Y))),
                 cbind(as.matrix(X), as.matrix(Y)))
    expect_equal(unname(as.matrix(cbind(X, as.csc.matrix(Y)))),
                 cbind(as.matrix(X), as.matrix(Y)))
})

test_that("cbind vectors", {
    set.seed(1)
    X <- as.coo.matrix(rsparsematrix(100, 50, .3))
    v <- as.sparse.vector(rnorm(100))
    dvec <- as.numeric(v)

    expect_s4_class(cbind(v, v), "dgCMatrix")
    expect_s4_class(cbind(as.csc.matrix(X), v), "dgCMatrix")
    expect_s4_class(cbind(v, as.csc.matrix(X)), "dgCMatrix")
    expect_s4_class(cbind(X, v), "dgTMatrix")
    expect_s4_class(cbind(v, X), "dgTMatrix")

    expect_equal(unname(as.matrix(cbind(v, v))),
                 unname(cbind(dvec, dvec)))
    expect_equal(unname(as.matrix(cbind(as.csc.matrix(X), v))),
                 unname(cbind(as.matrix(X), dvec)))
    expect_equal(unname(as.matrix(cbind(v, as.csc.matrix(X)))),
                 unname(cbind(dvec, as.matrix(X))))
    expect_equal(unname(as.matrix(cbind(X, v))),
                 unname(cbind(as.matrix(X), dvec)))
    expect_equal(unname(as.matrix(cbind(v, X))),
                 unname(cbind(dvec, as.matrix(X))))
})
