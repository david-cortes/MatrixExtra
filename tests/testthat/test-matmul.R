library("testthat")
library("Matrix")
library("MatrixExtra")
library("float")
context("Matrix multiplications")

set_new_matrix_behavior()
options("MatrixExtra.nthreads" = 2)

### TODO: add tests about the names of objects

test_that("matmult dense-CSC", {
    set.seed(1)
    A <- rsparsematrix(100, 50, .4)
    B <- rsparsematrix(50, 20, .4)
    A1 <- rsparsematrix(1, 50, .4)
    B1 <- rsparsematrix(50, 1, .4)

    expect_equal(as.matrix(A) %*% as.csc.matrix(B), as.matrix(A) %*% as.matrix(B))
    expect_equal(float::dbl(float::fl(as.matrix(A)) %*% as.csc.matrix(B)),
                 as.matrix(A) %*% as.matrix(B), tolerance=1e-5)
    expect_s4_class(float::fl(as.matrix(A)) %*% as.csc.matrix(B), "float32")

    expect_equal(as.matrix(A1) %*% as.csc.matrix(B), as.matrix(A1) %*% as.matrix(B))
    expect_equal(as.matrix(A) %*% as.csc.matrix(B1), as.matrix(A) %*% as.matrix(B1))
    expect_equal(as.matrix(A1) %*% as.csc.matrix(B1), as.matrix(A1) %*% as.matrix(B1))

    expect_equal(as.matrix(A) %*% as.csc.matrix(B, binary=TRUE),
                 as.matrix(A) %*% as.matrix(as.csc.matrix(B, binary=TRUE)))
    expect_equal(as.matrix(A) %*% as.csc.matrix(B, logical=TRUE),
                 as.matrix(A) %*% as.matrix(as.csc.matrix(B, logical=TRUE)))

    sy <- sparseMatrix(i= c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7),
                       symmetric=TRUE, dimnames = list(NULL, letters[1:7]))
    expect_s4_class(sy, "dsCMatrix")
    set.seed(1)
    sy_counterpart <- rsparsematrix(20, nrow(sy), .4)
    sy_counterpart <- as.matrix(sy_counterpart)

    expect_equal(sy_counterpart %*% sy, sy_counterpart %*% as.matrix(sy))

    tri <- matrix(c(1,2,0,4, 0,0,6,7, 0,0,8,9, 0,0,0,0), byrow=TRUE, nrow=4)
    tri <- as(tri, "triangularMatrix")
    expect_s4_class(tri, "dtCMatrix")
    set.seed(1)
    tri_counterpart <- rsparsematrix(20, nrow(tri), .4)
    tri_counterpart <- as.matrix(tri_counterpart)

    expect_equal(tri_counterpart %*% tri, tri_counterpart %*% as.matrix(tri))
})

test_that("tcrossprod dense-CSR", {
    set.seed(1)
    A <- rsparsematrix(100, 50, .4)
    B <- rsparsematrix(20, 50, .4)
    A1 <- rsparsematrix(1, 50, .4)
    B1 <- rsparsematrix(1, 50, .4)

    expect_equal(tcrossprod(as.matrix(A), as.csr.matrix(B)),
                 tcrossprod(as.matrix(A), as.matrix(as.csr.matrix(B))))
    expect_equal(float::dbl(tcrossprod(float::fl(as.matrix(A)), as.csr.matrix(B))),
                 tcrossprod(as.matrix(A), as.matrix(as.csr.matrix(B))),
                 tolerance=1e-5)
    expect_s4_class(tcrossprod(float::fl(as.matrix(A)), as.csr.matrix(B)), "float32")


    expect_equal(tcrossprod(as.matrix(A1), as.csr.matrix(B)),
                 tcrossprod(as.matrix(A1), as.matrix(as.csr.matrix(B))))
    expect_equal(tcrossprod(as.matrix(A), as.csr.matrix(B1)),
                 tcrossprod(as.matrix(A), as.matrix(as.csr.matrix(B1))))
    expect_equal(tcrossprod(as.matrix(A1), as.csr.matrix(B1)),
                 tcrossprod(as.matrix(A1), as.matrix(as.csr.matrix(B1))))

    sy <- sparseMatrix(i= c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7),
                       symmetric=TRUE, dimnames = list(NULL, letters[1:7]))
    sy <- as(sy, "RsparseMatrix")
    expect_s4_class(sy, "dsRMatrix")
    set.seed(1)
    sy_counterpart <- rsparsematrix(20, nrow(sy), .4)
    sy_counterpart <- as.matrix(sy_counterpart)

    expect_equal(tcrossprod(sy_counterpart, sy),
                 tcrossprod(sy_counterpart, as.matrix(sy)))

    tri <- matrix(c(1,2,0,4, 0,0,6,7, 0,0,8,9, 0,0,0,0), byrow=TRUE, nrow=4)
    tri <- as(tri, "triangularMatrix")
    tri <- as(tri, "RsparseMatrix")
    expect_s4_class(tri, "dtRMatrix")
    set.seed(1)
    tri_counterpart <- rsparsematrix(20, nrow(tri), .4)
    tri_counterpart <- as.matrix(tri_counterpart)

    expect_equal(tcrossprod(tri_counterpart, tri),
                 tcrossprod(tri_counterpart, as.matrix(tri)))
})

test_that("crossprod dense-CSC", {
    set.seed(1)
    A <- rsparsematrix(50, 100, .4)
    B <- rsparsematrix(50, 20, .4)

    expect_equal(crossprod(as.matrix(A), as.csc.matrix(B)),
                 crossprod(as.matrix(A), as.matrix(B)),
                 tolerance=1e-3)
})

test_that("matmult CSR-dense", {
    set.seed(1)
    A <- rsparsematrix(100, 50, .4)
    B <- rsparsematrix(50, 20, .4)

    expect_equal(as.csr.matrix(A) %*% as.matrix(B), as.matrix(A) %*% as.matrix(B))
})

test_that("tcrossprod CSR-dense", {
    set.seed(1)
    A <- rsparsematrix(100, 50, .4)
    B <- rsparsematrix(20, 50, .4)

    expect_equal(tcrossprod(as.csr.matrix(A), as.matrix(B)),
                 tcrossprod(as.matrix(A), as.matrix(B)))
})

test_that("tcrossprod dense-CSR", {
    set.seed(1)
    A <- rsparsematrix(100, 50, .4)
    B <- rsparsematrix(20, 50, .4)

    expect_equal(tcrossprod(as.matrix(A), as.csr.matrix(B)),
                 tcrossprod(as.matrix(A), as.matrix(B)))
})

test_that("matmult CSR-vector", {
    set.seed(1)
    A <- rsparsematrix(100, 50, .4)
    v <- rnorm(50)

    dvec <- as(v, "dsparseVector")
    ivec <- as(dvec, "isparseVector")
    lvec <- as(dvec, "lsparseVector")
    nvec <- as(dvec, "nsparseVector")

    num <- as.numeric(v)
    int <- as.integer(num)
    bool <- as.logical(num)

    lst_inputs <- list(
        dvec, ivec, lvec, nvec,
        num, int, bool
    )
    for (inp in lst_inputs) {
        expect_equal(as.csr.matrix(A) %*% inp,
                     as.matrix(A) %*% as.numeric(inp))
        expect_equal(unname(as.matrix(as.csr.matrix(A[,1,drop=FALSE]) %*% inp)),
                     unname(as.matrix(A[,1,drop=FALSE]) %*% as.numeric(inp)))
        expect_equal(unname(as.matrix(as.csr.matrix(A[1:10,1,drop=FALSE]) %*% inp)),
                     unname(as.matrix(A[1:10,1,drop=FALSE]) %*% as.numeric(inp)))
    }

    v[4] <- NA_real_
    A <- as.csr.matrix(A)
    A@x[10] <- NA_real_
    expect_equal(A %*% v, unname(as.matrix(as.csc.matrix(A) %*% v)))
})

test_that("float32 vectors", {
    set.seed(1)
    A <- rsparsematrix(100, 50, .4)
    v <- float::fl(rnorm(100))
    v2 <- float::fl(rnorm(50))

    expect_equal(float::dbl(v %*% A),
                 float::dbl(v) %*% as.matrix(A),
                 tolerance=1e-5)

    ### sparse CSC
    expect_equal(unname(as.matrix(v2 %*% as.csc.matrix(A[1,,drop=FALSE]))),
                 float::dbl(v2) %*% as.matrix(A)[1,,drop=FALSE],
                 tolerance=1e-5)

    ### sparse CSC
    expect_equal(unname(as.matrix(v2 %*% as.csc.matrix(A[1,1:10,drop=FALSE]))),
                 float::dbl(v2) %*% as.matrix(A)[1,1:10,drop=FALSE],
                 tolerance=1e-5)

    expect_equal(float::dbl(tcrossprod(v2, as.csr.matrix(A))),
                 tcrossprod(float::dbl(v2), as.matrix(A)),
                 tolerance=1e-5)

    expect_equal(float::dbl(tcrossprod(v2, as.csr.matrix(A[1:10,]))),
                 tcrossprod(float::dbl(v2), as.matrix(A)[1:10,]),
                 tolerance=1e-5)

    ### sparse CSC
    expect_equal(unname(as.matrix(tcrossprod(v, as.csr.matrix(A[,1,drop=FALSE])))),
                 tcrossprod(float::dbl(v), as.matrix(A[,1,drop=FALSE])),
                 tolerance=1e-5)

    ### sparse CSC
    expect_equal(unname(as.matrix(tcrossprod(v, as.csr.matrix(A[1:10,1,drop=FALSE])))),
                 tcrossprod(float::dbl(v), as.matrix(A[1:10,1,drop=FALSE])),
                 tolerance=1e-5)

    expect_equal(float::dbl(crossprod(v, as.csc.matrix(A))),
                 crossprod(float::dbl(v), as.matrix(A)),
                 tolerance=1e-5)

    expect_equal(float::dbl(crossprod(v, as.csc.matrix(A[,1:10,drop=FALSE]))),
                 crossprod(float::dbl(v), as.matrix(A)[,1:10,drop=FALSE]),
                 tolerance=1e-5)

    expect_equal(float::dbl(as.csr.matrix(A) %*% v2),
                 as.matrix(A) %*% float::dbl(v2),
                 tolerance=1e-5)

    # sparse CSR
    expect_equal(unname(as.matrix(as.csr.matrix(A[,1,drop=FALSE]) %*% v2)),
                 as.matrix(A[,1,drop=FALSE]) %*% float::dbl(v2),
                 tolerance=1e-5)

    # sparse CSR
    expect_equal(unname(as.matrix(as.csr.matrix(A[1:10,1,drop=FALSE]) %*% v2)),
                 as.matrix(A[1:10,1,drop=FALSE]) %*% float::dbl(v2),
                 tolerance=1e-5)
})

set_new_matrix_behavior()
