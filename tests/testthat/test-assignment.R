library("testthat")
library("Matrix")
library("MatrixExtra")
context("Assignment")


X1d <- matrix(0, nrow=4, ncol=3)
X2d <- matrix(1, nrow=4, ncol=3)
set.seed(1)
X3d <- matrix(rnorm(12), nrow=4, ncol=3)
X4d <- matrix(c(0,0,0, 1,2,3, 0,0,0, 0,1,0), nrow=4, ncol=3, byrow=TRUE)

X1 <- as.csr.matrix(X1d)
X2 <- as.csr.matrix(X2d)
X3 <- as.csr.matrix(X3d)
X4 <- as.csr.matrix(X4d)

set.seed(1)
X5 <- rsparsematrix(111, 23, .1, repr="R")
X6 <- rsparsematrix(222, 6, .1, repr="R")
X5d <- as.matrix(X5)
X6d <- as.matrix(X6)

get_sparse_vector <- function(size) {
    set.seed(1)
    v <- rsparsematrix(size, 1, .75, repr="T")
    v <- as.sparse.vector(v)
    return(v)
}

test_that("Set single to zero", {

    run_tests <- function(X, Xd) {
        X_orig <- X
        Xd_orig <- Xd

        X[2,] <- 0
        Xd[2,] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[,2] <- 0
        Xd[,2] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[2,3] <- 0
        Xd[2,3] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[3,2] <- 0
        Xd[3,2] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[4,1] <- 0
        Xd[4,1] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[4,2] <- 0
        Xd[4,2] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set single to constant", {
    run_tests <- function(X, Xd) {
        X_orig <- X
        Xd_orig <- Xd

        X[2,] <- 111
        Xd[2,] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[,2] <- 111
        Xd[,2] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[2,3] <- 111
        Xd[2,3] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[3,2] <- 111
        Xd[3,2] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[4,1] <- 111
        Xd[4,1] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[4,2] <- 111
        Xd[4,2] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set single to dense vector", {
    run_tests <- function(X, Xd) {
        X_orig <- X
        Xd_orig <- Xd

        X[2,] <- seq(1, ncol(X))
        Xd[2,] <- seq(1, ncol(Xd))
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[3,] <- seq(1, ncol(X))
        Xd[3,] <- seq(1, ncol(Xd))
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[,2] <- seq(1, nrow(X))
        Xd[,2] <- seq(1, nrow(Xd))
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[,3] <- seq(1, nrow(X))
        Xd[,3] <- seq(1, nrow(Xd))
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        row_div2 <- nrow(X) > 2 && (nrow(X) %% 2) == 0
        row_div3 <- nrow(X) > 3 && (nrow(X) %% 3) == 0
        col_div2 <- ncol(X) > 2 && (ncol(X) %% 2) == 0
        col_div3 <- ncol(X) > 3 && (ncol(X) %% 3) == 0

        if (col_div2) {
            X[2,] <- seq(1, ncol(X)/2)
            Xd[2,] <- seq(1, ncol(Xd)/2)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig

            X[3,] <- seq(1, ncol(X)/2)
            Xd[3,] <- seq(1, ncol(Xd)/2)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }

        if (col_div3) {
            X[2,] <- seq(1, ncol(X)/3)
            Xd[2,] <- seq(1, ncol(Xd)/3)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig

            X[3,] <- seq(1, ncol(X)/3)
            Xd[3,] <- seq(1, ncol(Xd)/3)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }

        if (row_div2) {
            X[,2] <- seq(1, nrow(X)/2)
            Xd[,2] <- seq(1, nrow(Xd)/2)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig

            X[,3] <- seq(1, nrow(X)/2)
            Xd[,3] <- seq(1, nrow(Xd)/2)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }

        if (row_div3) {
            X[,2] <- seq(1, nrow(X)/3)
            Xd[,2] <- seq(1, nrow(Xd)/3)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig

            X[,3] <- seq(1, nrow(X)/3)
            Xd[,3] <- seq(1, nrow(Xd)/3)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set single to sparse vector", {
    assign_bytype_row <- function(X, Xd, X_orig, Xd_orig, v, i) {
        X[i,] <- as.sparse.vector(v)
        Xd[i,] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[i,] <- as.csr.matrix(v)
        Xd[i,] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[i,] <- as.csc.matrix(v)
        Xd[i,] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[i,] <- as.coo.matrix(v)
        Xd[i,] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
    }

    assign_bytype_col <- function(X, Xd, X_orig, Xd_orig, v, j) {
        X[,j] <- as.sparse.vector(v)
        Xd[,j] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[,j] <- as.csr.matrix(v)
        Xd[,j] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[,j] <- as.csc.matrix(v)
        Xd[,j] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[,j] <- as.coo.matrix(v)
        Xd[,j] <- as.numeric(v)
        expect_equal(unname(as.matrix(X)), unname(Xd))
    }


    run_tests <- function(X, Xd) {
        X_orig <- X
        Xd_orig <- Xd

        v <- get_sparse_vector(ncol(X))
        assign_bytype_row(X, Xd, X_orig, Xd_orig, v, 2)

        v <- get_sparse_vector(ncol(X))
        assign_bytype_row(X, Xd, X_orig, Xd_orig, v, 3)

        v <- get_sparse_vector(nrow(X))
        assign_bytype_col(X, Xd, X_orig, Xd_orig, v, 2)

        v <- get_sparse_vector(nrow(X))
        assign_bytype_col(X, Xd, X_orig, Xd_orig, v, 3)

        row_div2 <- nrow(X) > 2 && (nrow(X) %% 2) == 0
        row_div3 <- nrow(X) > 3 && (nrow(X) %% 3) == 0
        col_div2 <- ncol(X) > 2 && (ncol(X) %% 2) == 0
        col_div3 <- ncol(X) > 3 && (ncol(X) %% 3) == 0

        if (col_div2) {
            v <- get_sparse_vector(ncol(X)/2)
            assign_bytype_row(X, Xd, X_orig, Xd_orig, v, 2)

            v <- get_sparse_vector(ncol(X)/2)
            assign_bytype_row(X, Xd, X_orig, Xd_orig, v, 3)
        }

        if (col_div3) {
            v <- get_sparse_vector(ncol(X)/3)
            assign_bytype_row(X, Xd, X_orig, Xd_orig, v, 2)

            v <- get_sparse_vector(ncol(X)/3)
            assign_bytype_row(X, Xd, X_orig, Xd_orig, v, 3)
        }

        if (row_div2) {
            v <- get_sparse_vector(nrow(X)/2)
            assign_bytype_col(X, Xd, X_orig, Xd_orig, v, 2)

            v <- get_sparse_vector(nrow(X)/2)
            assign_bytype_col(X, Xd, X_orig, Xd_orig, v, 3)
        }

        if (row_div3) {
            v <- get_sparse_vector(nrow(X)/3)
            assign_bytype_col(X, Xd, X_orig, Xd_orig, v, 2)

            v <- get_sparse_vector(nrow(X)/3)
            assign_bytype_col(X, Xd, X_orig, Xd_orig, v, 3)
        }
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set row sequence to constant", {
    s1 <- seq(1L, 2L)
    s2 <- seq(2L, 4L)
    s3 <- seq(3L, 1L)

    run_tests <- function(X, Xd, s) {
        X_orig <- X
        Xd_orig <- Xd

        s_large <- seq(2L, nrow(X) - 1L)

        X[s1, ] <- 0
        Xd[s1, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s2, ] <- 0
        Xd[s2, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s3, ] <- 0
        Xd[s3, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s_large, ] <- 0
        Xd[s_large, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s1, ] <- 111
        Xd[s1, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s2, ] <- 111
        Xd[s2, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s3, ] <- 111
        Xd[s3, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s_large, ] <- 111
        Xd[s_large, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set column sequence to constant", {
    s1 <- seq(1L, 2L)
    s2 <- seq(2L, 3L)

    run_tests <- function(X, Xd) {
        X_orig <- X
        Xd_orig <- Xd

        s_large <- seq(2L, ncol(X) - 1L)

        X[, s1] <- 0
        Xd[, s1] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s2] <- 0
        Xd[, s2] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s_large] <- 0
        Xd[, s_large] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s1] <- 111
        Xd[, s1] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s2] <- 111
        Xd[, s2] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s_large] <- 111
        Xd[, s_large] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set random rows to constant", {

    run_tests <- function(X, Xd) {
        s1 <- c(1L, nrow(X))
        s2 <- c(1L, 4L, 3L)
        s3 <- c(TRUE, FALSE)
        s4 <- c(FALSE, TRUE)
        set.seed(1)
        s5 <- sample(nrow(X), floor(nrow(X)/2), replace=FALSE)

        X_orig <- X
        Xd_orig <- Xd

        X[s1, ] <- 0
        Xd[s1, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s2, ] <- 0
        Xd[s2, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s3, ] <- 0
        Xd[s3, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s4, ] <- 0
        Xd[s4, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s5, ] <- 0
        Xd[s5, ] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s1, ] <- 111
        Xd[s1, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s2, ] <- 111
        Xd[s2, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s3, ] <- 111
        Xd[s3, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s4, ] <- 111
        Xd[s4, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[s5, ] <- 111
        Xd[s5, ] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set random columns to constant", {
    run_tests <- function(X, Xd) {
        s1 <- c(1L, ncol(X))
        s2 <- c(1L, 3L)
        s3 <- c(TRUE, FALSE)
        s4 <- c(FALSE, TRUE)
        set.seed(1)
        s5 <- sample(ncol(X), floor(ncol(X)/2), replace=FALSE)

        X_orig <- X
        Xd_orig <- Xd

        X[, s1] <- 0
        Xd[, s1] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s2] <- 0
        Xd[, s2] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s3] <- 0
        Xd[, s3] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s4] <- 0
        Xd[, s4] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s5] <- 0
        Xd[, s5] <- 0
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s1] <- 111
        Xd[, s1] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s2] <- 111
        Xd[, s2] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s3] <- 111
        Xd[, s3] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s4] <- 111
        Xd[, s4] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig

        X[, s5] <- 111
        Xd[, s5] <- 111
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set single col in random rows to const", {
    run_tests <- function(X, Xd) {
        s1 <- c(1L, nrow(X))
        s2 <- c(1L, 3L)
        s3 <- c(TRUE, FALSE)
        s4 <- c(FALSE, TRUE)
        set.seed(1)
        s5 <- sample(nrow(X), floor(nrow(X)/2), replace=FALSE)

        X_orig <- X
        Xd_orig <- Xd

        stry <- list(s1, s2, s3, s4, s5)
        jtry <- seq(1L, ncol(X))
        for (s in stry) {
            for (j in jtry) {
                X[s, j] <- 0
                Xd[s, j] <- 0
                expect_equal(unname(as.matrix(X)), unname(Xd))
                X <- X_orig
                Xd <- Xd_orig

                X[s, j] <- 111
                Xd[s, j] <- 111
                expect_equal(unname(as.matrix(X)), unname(Xd))
                X <- X_orig
                Xd <- Xd_orig
            }
        }
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set random cols in single row to const", {
    run_tests <- function(X, Xd) {
        s1 <- c(1L, ncol(X))
        s2 <- c(1L, 3L)
        s3 <- c(TRUE, FALSE)
        s4 <- c(FALSE, TRUE)
        set.seed(1)
        s5 <- sample(ncol(X), floor(ncol(X)/2), replace=FALSE)

        X_orig <- X
        Xd_orig <- Xd

        stry <- list(s1, s2, s3, s4, s5)
        itry <- seq(1L, ncol(X))
        for (s in stry) {
            for (i in itry) {
                X[i, s] <- 0
                Xd[i, s] <- 0
                expect_equal(unname(as.matrix(X)), unname(Xd))
                X <- X_orig
                Xd <- Xd_orig

                X[i, s] <- 111
                Xd[i, s] <- 111
                expect_equal(unname(as.matrix(X)), unname(Xd))
                X <- X_orig
                Xd <- Xd_orig
            }
        }
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set arbitrary rows and columns to constant", {
    run_tests <- function(X, Xd) {
        i1 <- c(1L, nrow(X))
        i2 <- c(1L, 3L)
        i3 <- c(TRUE, FALSE)
        i4 <- c(FALSE, TRUE)
        set.seed(1)
        i5 <- sample(nrow(X), floor(nrow(X)/2), replace=FALSE)

        j1 <- c(1L, ncol(X))
        j2 <- c(1L, 3L)
        j3 <- c(TRUE, FALSE)
        j4 <- c(FALSE, TRUE)
        set.seed(1)
        j5 <- sample(ncol(X), floor(ncol(X)/2), replace=FALSE)

        X_orig <- X
        Xd_orig <- Xd

        itry <- list(i1, i2, i3, i4, i5)
        jtry <- list(j1, j2, j3, j4, j5)
        for (i in itry) {
            for (j in jtry) {
                X[i, j] <- 0
                Xd[i, j] <- 0
                expect_equal(unname(as.matrix(X)), unname(Xd))
                X <- X_orig
                Xd <- Xd_orig

                X[i, j] <- 111
                Xd[i, j] <- 111
                expect_equal(unname(as.matrix(X)), unname(Xd))
                X <- X_orig
                Xd <- Xd_orig
            }
        }
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set row sequence to sparse matrix", {
    run_tests <- function(X, Xd) {
        s1 <- 1:2
        s2 <- 2:3
        s3 <- 2:1
        s4 <- 3:2
        s5 <- seq(1, nrow(X))
        s6 <- seq(2, nrow(X)-1)
        s7 <- seq(3, nrow(X))

        X_orig <- X
        Xd_orig <- Xd
        s_ <- list(s1, s2, s3, s4, s5, s6, s7)
        s_ <- s_[sapply(s_, function(x) length(x) > 1)]
        for (s in s_) {
            set.seed(111)
            Y <- rsparsematrix(length(s), ncol(X), .5)

            X[s, ] <- as.csr.matrix(Y)
            Xd[s, ] <- as.matrix(Y)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig

            X[s, ] <- as.csc.matrix(Y)
            Xd[s, ] <- as.matrix(Y)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig

            X[s, ] <- as.coo.matrix(Y)
            Xd[s, ] <- as.matrix(Y)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }
    }

    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set random rows to sparse matrix", {
    
    run_tests <- function(X, Xd) {
        X_orig <- X
        Xd_orig <- Xd
        
        assign_values <- function(X, Xd, Y, s, X_orig, Xd_orig) {
            X[s, ] <- as.coo.matrix(Y)
            Xd[s, ] <- as.matrix(Y)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
            
            X[s, ] <- as.csr.matrix(Y)
            Xd[s, ] <- as.matrix(Y)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
            
            X[s, ] <- as.csc.matrix(Y)
            Xd[s, ] <- as.matrix(Y)
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }
        
        set.seed(111)
        s <- sample(nrow(X), nrow(X) / 2, replace=FALSE)
        s <- s[order(s)]
        Y <- rsparsematrix(length(s), ncol(X), .5)
        assign_values(X, Xd, Y, s, X_orig, Xd_orig)
        
        set.seed(112)
        s <- sample(nrow(X), nrow(X), replace=FALSE)
        s <- s[order(s)]
        Y <- rsparsematrix(length(s), ncol(X), .5)
        assign_values(X, Xd, Y, s, X_orig, Xd_orig)
        
        set.seed(113)
        s <- sample(nrow(X), 2, replace=FALSE)
        s <- s[order(s)]
        Y <- rsparsematrix(length(s), ncol(X), .5)
        assign_values(X, Xd, Y, s, X_orig, Xd_orig)
    }
    
    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})

test_that("Set full sequences to zero", {
    run_tests <- function(X, Xd) {
        v1 <- 0
        v2 <- numeric(2)
        v3 <- as.sparse.vector(emptySparse(nrow=1, ncol=2))
        v4 <- emptySparse(nrow=1, ncol=3)
        
        v3d <- as.numeric(v3)
        v4d <- as.matrix(v4)
        
        X_orig <- X
        Xd_orig <- Xd
        
        suppressWarnings(X[1:nrow(X), 1:ncol(X)] <- v1)
        Xd[1:nrow(Xd), 1:ncol(Xd)] <- v1
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
        
        suppressWarnings(X[rev(1:nrow(X)), rev(1:ncol(X))] <- v1)
        Xd[rev(1:nrow(Xd)), rev(1:ncol(Xd))] <- v1
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
        
        if ((nrow(X) %% 2) == 0 || (ncol(X) %% 2) == 0) {
            suppressWarnings(X[1:nrow(X), 1:ncol(X)] <- v2)
            Xd[1:nrow(Xd), 1:ncol(Xd)] <- v2
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
            
            suppressWarnings(X[rev(1:nrow(X)), rev(1:ncol(X))] <- v2)
            Xd[rev(1:nrow(Xd)), rev(1:ncol(Xd))] <- v2
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
            
            suppressWarnings(X[1:nrow(X), 1:ncol(X)] <- v3)
            Xd[1:nrow(Xd), 1:ncol(Xd)] <- v3d
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
            
            suppressWarnings(X[rev(1:nrow(X)), rev(1:ncol(X))] <- v3)
            Xd[rev(1:nrow(Xd)), rev(1:ncol(Xd))] <- v3d
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }
        
        if (ncol(X) > 6) {
            X[, 1:6] <- v2
            Xd[, 1:6] <- v2
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
            
            X[, 1:6] <- v3
            Xd[, 1:6] <- v3d
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
            
            X[, 1:6] <- v4
            Xd[, 1:6] <- v4d
            expect_equal(unname(as.matrix(X)), unname(Xd))
            X <- X_orig
            Xd <- Xd_orig
        }
        
        expect_error(X[2, 1:ncol(X)] <- emptySparse(nrow=3, ncol=ncol(X)))
        
        X[2, 1:ncol(X)] <- emptySparse(nrow=1, ncol=ncol(X))
        Xd[2, 1:ncol(Xd)] <- as.matrix(emptySparse(nrow=1, ncol=ncol(Xd)))
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
        
        X[1:3, 1:ncol(X)] <- emptySparse(nrow=3, ncol=ncol(X))
        Xd[1:3, 1:ncol(Xd)] <- as.matrix(emptySparse(nrow=3, ncol=ncol(Xd)))
        expect_equal(unname(as.matrix(X)), unname(Xd))
        X <- X_orig
        Xd <- Xd_orig
    }
    
    run_tests(X1, X1d)
    run_tests(X2, X2d)
    run_tests(X3, X3d)
    run_tests(X4, X4d)
    run_tests(X5, X5d)
    run_tests(X6, X6d)
})
