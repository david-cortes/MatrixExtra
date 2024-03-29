library("testthat")
library("Matrix")
library("MatrixExtra")
context("Rbinding matrices")

test_that("Types of resulting objects", {
    set.seed(1)
    Mat <- rsparsematrix(10, 3, .4)
    Mat <- as.csr.matrix(Mat)
    vec <- as(Mat[3, , drop=FALSE], "sparseVector")

    expect_s4_class(rbind2(Mat, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, Mat, Mat), "dgRMatrix")
    expect_s4_class(rbind2(vec, vec), "dgRMatrix")

    expect_s4_class(rbind2(Mat, vec), "dgRMatrix")
    expect_s4_class(rbind(Mat, vec), "dgRMatrix")
    expect_s4_class(rbind2(vec, Mat), "dgRMatrix")
    expect_s4_class(rbind(vec, Mat), "dgRMatrix")
    expect_s4_class(rbind(vec, Mat, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, vec, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, Mat, vec), "dgRMatrix")
    
    dvec <- as.numeric(vec)
    ivec <- as.integer(dvec)
    lvec <- as.logical(dvec)
    
    expect_s4_class(rbind2(Mat, dvec), "dgRMatrix")
    expect_s4_class(rbind(Mat, dvec), "dgRMatrix")
    expect_s4_class(rbind2(dvec, Mat), "dgRMatrix")
    expect_s4_class(rbind(dvec, Mat), "dgRMatrix")
    expect_s4_class(rbind(dvec, Mat, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, dvec, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, Mat, dvec), "dgRMatrix")
    
    expect_s4_class(rbind2(Mat, ivec), "dgRMatrix")
    expect_s4_class(rbind(Mat, ivec), "dgRMatrix")
    expect_s4_class(rbind2(ivec, Mat), "dgRMatrix")
    expect_s4_class(rbind(ivec, Mat), "dgRMatrix")
    expect_s4_class(rbind(ivec, Mat, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, dvec, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, Mat, ivec), "dgRMatrix")
    
    expect_s4_class(rbind2(Mat, lvec), "dgRMatrix")
    expect_s4_class(rbind(Mat, lvec), "dgRMatrix")
    expect_s4_class(rbind2(lvec, Mat), "dgRMatrix")
    expect_s4_class(rbind(lvec, Mat), "dgRMatrix")
    expect_s4_class(rbind(lvec, Mat, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, lvec, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, Mat, lvec), "dgRMatrix")

    MatBin <- as.csr.matrix(Mat, binary=TRUE)
    vecBin <- as(vec, "nsparseVector")

    expect_s4_class(rbind2(MatBin, MatBin), "ngRMatrix")
    expect_s4_class(rbind(MatBin, MatBin, MatBin), "ngRMatrix")
    expect_s4_class(rbind2(vecBin, vecBin), "ngRMatrix")

    expect_s4_class(rbind2(MatBin, vecBin), "ngRMatrix")
    expect_s4_class(rbind(MatBin, vecBin), "ngRMatrix")
    expect_s4_class(rbind2(vecBin, MatBin), "ngRMatrix")
    expect_s4_class(rbind(vecBin, MatBin), "ngRMatrix")
    expect_s4_class(rbind(vecBin, MatBin, MatBin), "ngRMatrix")
    expect_s4_class(rbind(MatBin, vecBin, MatBin), "ngRMatrix")
    expect_s4_class(rbind(MatBin, MatBin, vecBin), "ngRMatrix")

    expect_s4_class(rbind(Mat, MatBin), "dgRMatrix")
    expect_s4_class(rbind(MatBin, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, vecBin), "dgRMatrix")
    expect_s4_class(rbind(MatBin, vec), "dgRMatrix")
    expect_s4_class(rbind(vec, MatBin), "dgRMatrix")
    expect_s4_class(rbind(vecBin, Mat), "dgRMatrix")
    expect_s4_class(rbind(vecBin, vec), "dgRMatrix")
    expect_s4_class(rbind(vec, vecBin), "dgRMatrix")

    MatBool <- as.csr.matrix(Mat, logical=TRUE)
    vecBool <- as(vec, "lsparseVector")

    expect_s4_class(rbind2(MatBool, MatBool), "lgRMatrix")
    expect_s4_class(rbind(MatBool, MatBool, MatBool), "lgRMatrix")
    expect_s4_class(rbind2(vecBool, vecBool), "lgRMatrix")

    expect_s4_class(rbind2(MatBool, vecBool), "lgRMatrix")
    expect_s4_class(rbind(MatBool, vecBool), "lgRMatrix")
    expect_s4_class(rbind2(vecBool, MatBool), "lgRMatrix")
    expect_s4_class(rbind(vecBool, MatBool), "lgRMatrix")
    expect_s4_class(rbind(vecBool, MatBool, MatBool), "lgRMatrix")
    expect_s4_class(rbind(MatBool, vecBool, MatBool), "lgRMatrix")
    expect_s4_class(rbind(MatBool, MatBool, vecBool), "lgRMatrix")

    expect_s4_class(rbind(Mat, MatBool), "dgRMatrix")
    expect_s4_class(rbind(MatBool, Mat), "dgRMatrix")
    expect_s4_class(rbind(Mat, vecBool), "dgRMatrix")
    expect_s4_class(rbind(MatBool, vec), "dgRMatrix")
    expect_s4_class(rbind(vec, MatBool), "dgRMatrix")
    expect_s4_class(rbind(vecBool, Mat), "dgRMatrix")
    expect_s4_class(rbind(vecBool, vec), "dgRMatrix")
    expect_s4_class(rbind(vec, vecBool), "dgRMatrix")
})

test_that("Contents of resulting object", {
    set.seed(1)
    X <- rsparsematrix(5, 3, .4)
    X <- as.csr.matrix(X)
    Y <- rsparsematrix(2, 3, .4)
    Y <- as.csr.matrix(Y)
    v <- as(X[2, , drop=FALSE], "sparseVector")

    expect_equal(unname(as.matrix(rbind(X, Y))), rbind(as.matrix(X), as.matrix(Y)))
    expect_equal(unname(as.matrix(rbind(X, v))), rbind(as.matrix(X), as.numeric(v)))
    expect_equal(unname(as.matrix(rbind(v, Y))), rbind(as.numeric(v), as.matrix(Y)))

    expect_equal(unname(as.matrix(rbind(X, as.numeric(v)))), rbind(as.matrix(X), as.numeric(v)))
    expect_equal(unname(as.matrix(rbind(as.numeric(v), Y))), rbind(as.numeric(v), as.matrix(Y)))

    X.mat <- as.matrix(X)
    X.mat[1,1] <- NA
    mode(X.mat) <- "logical"
    X.mat <- as(X.mat, "RsparseMatrix")

    expect_equal(unname(as.matrix(rbind(X.mat, Y))), rbind(as.matrix(X.mat), as.matrix(Y)))
    expect_equal(unname(as.matrix(rbind(X.mat, v))), rbind(as.matrix(X.mat), as.numeric(v)))
    expect_equal(unname(as.matrix(rbind(v, X.mat))), rbind(as.numeric(v), as.matrix(X.mat)))
})

test_that("Batched rbinding", {
    set.seed(1)
    X <- rsparsematrix(5, 3, .4)
    X <- as.csr.matrix(X)
    v <- as(X[2, , drop=FALSE], "sparseVector")
    v_csr <- as.csr.matrix(v)

    expect_equal(rbind_csr(v, v, v), rbind(v_csr, v_csr, v_csr))
    expect_equal(rbind_csr(X, X, X), rbind(X, X, X))

    expect_equal(rbind_csr(v, X, X), rbind(v_csr, X, X))
    expect_equal(rbind_csr(X, v, X), rbind(X, v_csr, X))
    expect_equal(rbind_csr(X, X, v), rbind(X, X, v_csr))
    expect_equal(rbind_csr(X, v, v), rbind(X, v_csr, v_csr))
    expect_equal(rbind_csr(v, v, X), rbind(v_csr, v_csr, X))

    v_bool <- as(v, "lsparseVector")
    v_bool_csr <- as.csr.matrix(v_bool, logical=TRUE)
    X_bool <- as.csr.matrix(X, logical=TRUE)
    expect_equal(rbind_csr(v_bool, v_bool, v_bool), rbind(v_bool_csr, v_bool_csr, v_bool_csr))
    expect_equal(rbind_csr(X_bool, v_bool, v_bool), rbind(X_bool, v_bool_csr, v_bool_csr))

    v_bin <- as(v, "nsparseVector")
    v_bin_csr <- as.csr.matrix(v_bin, binary=TRUE)
    X_bin <- as.csr.matrix(X, binary=TRUE)
    expect_equal(rbind_csr(v_bin, v_bin, v_bin), rbind(v_bin_csr, v_bin_csr, v_bin_csr))
    expect_equal(rbind_csr(X_bin, v_bin, v_bin), rbind(X_bin, v_bin_csr, v_bin_csr))

    v_int <- as(v, "isparseVector")
    v_int_csr <- as.csr.matrix(v_int)
    expect_equal(rbind_csr(v_int, v_int, v_int), rbind(v_int_csr, v_int_csr, v_int_csr))

    expect_equal(rbind_csr(v, v_bin, v_bool), rbind(v_csr, v_bin_csr, v_bin_csr))
    expect_equal(rbind_csr(X, X_bin, X_bool), rbind(X, X_bin, X_bool))

    expect_equal(rbind_csr(X, v_bin, X_bin, X_bool, v_int), rbind(X, v_bin_csr, X_bin, X_bool, v_int_csr))
    
    X_dempty <- matrix(numeric(), ncol=ncol(X))
    X_rempty <- as.csr.matrix(X_dempty)
    X_cempty <- as.csc.matrix(X_dempty)
    
    expect_equal(unname(as.matrix(rbind_csr(X_dempty, X_cempty))),
                 unname(rbind(X_dempty, X_dempty)))
    expect_equal(unname(as.matrix(rbind_csr(X_cempty, X_cempty))),
                 unname(rbind(X_dempty, X_dempty)))
})

test_that("rbind COO", {
    set.seed(1)
    X <- as.coo.matrix(rsparsematrix(10, 5, .4))
    Y <- as.csr.matrix(rsparsematrix(5,  5, .4))
    Z <- as.csc.matrix(rsparsematrix(5,  5, .4))
    v <- as(X[2, , drop=FALSE], "sparseVector")
    dvec <- as.numeric(v)
    ivec <- as.integer(dvec)
    lvec <- as.logical(dvec)
    
    expect_s4_class(rbind(X, X), "TsparseMatrix")
    expect_s4_class(rbind2(X, X), "TsparseMatrix")
    expect_s4_class(rbind(X, Y), "RsparseMatrix")
    expect_s4_class(rbind(Y, X), "RsparseMatrix")
    expect_s4_class(rbind(X, Z), "TsparseMatrix")
    expect_s4_class(rbind(Z, X), "TsparseMatrix")
    
    expect_s4_class(rbind(X, v), "TsparseMatrix")
    expect_s4_class(rbind(v, X), "TsparseMatrix")
    expect_s4_class(rbind(X, dvec), "TsparseMatrix")
    expect_s4_class(rbind(dvec, X), "TsparseMatrix")
    expect_s4_class(rbind(X, ivec), "TsparseMatrix")
    expect_s4_class(rbind(ivec, X), "TsparseMatrix")
    expect_s4_class(rbind(X, lvec), "TsparseMatrix")
    expect_s4_class(rbind(lvec, X), "TsparseMatrix")
    
    
    expect_equal(unname(as.matrix(rbind(X, X))), rbind(as.matrix(X), as.matrix(X)))
    expect_equal(unname(as.matrix(rbind2(X, X))), rbind(as.matrix(X), as.matrix(X)))
    expect_equal(unname(as.matrix(rbind(X, Y))), rbind(as.matrix(X), as.matrix(Y)))
    expect_equal(unname(as.matrix(rbind(Y, X))), rbind(as.matrix(Y), as.matrix(X)))
    expect_equal(unname(as.matrix(rbind(X, Z))), rbind(as.matrix(X), as.matrix(Z)))
    expect_equal(unname(as.matrix(rbind(Z, X))), rbind(as.matrix(Z), as.matrix(X)))
    
    expect_equal(unname(as.matrix(rbind(X, v))), unname(rbind(as.matrix(X), dvec)))
    expect_equal(unname(as.matrix(rbind(v, X))), unname(rbind(dvec, as.matrix(X))))
    expect_equal(unname(as.matrix(rbind(X, dvec))), unname(rbind(as.matrix(X), dvec)))
    expect_equal(unname(as.matrix(rbind(dvec, X))), unname(rbind(dvec, as.matrix(X))))
    expect_equal(unname(as.matrix(rbind(X, ivec))), unname(rbind(as.matrix(X), ivec)))
    expect_equal(unname(as.matrix(rbind(ivec, X))), unname(rbind(ivec, as.matrix(X))))
    expect_equal(unname(as.matrix(rbind(X, lvec))), unname(rbind(as.matrix(X), lvec)))
    expect_equal(unname(as.matrix(rbind(lvec, X))), unname(rbind(lvec, as.matrix(X))))
})
