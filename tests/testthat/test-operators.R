library("testthat")
library("Matrix")
library("MatrixExtra")
context("Adding and multiplying matrices")

set.seed(1)
Mat1 <- rsparsematrix(100, 35, .4)
Mat2 <- rsparsematrix(100, 35, .6)
csr1 <- as.csr.matrix(Mat1)
csc1 <- as.csc.matrix(Mat1)
mat1 <- as.matrix(Mat1)
csr2 <- as.csr.matrix(Mat2)
csc2 <- as.csc.matrix(Mat2)
mat2 <- as.matrix(Mat2)
emat <- new("dgRMatrix")
emat@Dim <- csr1@Dim
emat@Dimnames <- csr1@Dimnames
emat@p <- integer(length(csr1@p))
eden <- as.matrix(emat)

expect_unmodified <- function(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden) {
    expect_equal(as.matrix(csr1), mat1)
    expect_equal(as.matrix(csr2), mat2)
    expect_equal(as.matrix(csc1), mat1)
    expect_equal(as.matrix(csc2), mat2)
    expect_equal(as.matrix(emat), eden)
}

test_that("Operations CSR-CSR", {
    expect_equal(as.matrix(csr1 + csr2), as.matrix(mat1 + mat2))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 + csc2), as.matrix(mat1 + mat2))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 + csr1), as.matrix(mat1 + mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 + emat), as.matrix(mat1 + eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat + emat), as.matrix(eden + eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 + csr2, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 + csc2, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 + csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 + emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat + emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)

    expect_equal(as.matrix(csr2 + csr1), as.matrix(mat2 + mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csc2 + csr1), as.matrix(mat2 + mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 + csr1), as.matrix(mat1 + mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat + csr1), as.matrix(eden + mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat + emat), as.matrix(eden + eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr2 + csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csc2 + csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 + csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat + csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat + emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)

    expect_equal(as.matrix(csr1 - csr2), as.matrix(mat1 - mat2))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 - csc2), as.matrix(mat1 - mat2))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 - csr1), as.matrix(mat1 - mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 - emat), as.matrix(mat1 - eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat - emat), as.matrix(eden - eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 - csr2, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 - csc2, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 - csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 - emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat - emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)

    expect_equal(as.matrix(csr2 - csr1), as.matrix(mat2 - mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csc2 - csr1), as.matrix(mat2 - mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 - csr1), as.matrix(mat1 - mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat - csr1), as.matrix(eden - mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat - emat), as.matrix(eden - eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr2 - csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csc2 - csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 - csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat - csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat - emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)

    expect_equal(as.matrix(csr1 * csr2), as.matrix(mat1 * mat2))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 * csc2), as.matrix(mat1 * mat2))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 * csr1), as.matrix(mat1 * mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 * emat), as.matrix(mat1 * eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat * emat), as.matrix(eden * eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 * csr2, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 * csc2, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 * csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 * emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat * emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)

    expect_equal(as.matrix(csr2 * csr1), as.matrix(mat2 * mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csc2 * csr1), as.matrix(mat2 * mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 * csr1), as.matrix(mat1 * mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat * csr1), as.matrix(eden * mat1))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat * emat), as.matrix(eden * eden))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr2 * csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csc2 * csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 * csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat * csr1, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat * emat, "dgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    
    expect_equal(as.matrix(csr1 | csr2), unname(as.matrix(mat1 | mat2)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 | csc2), unname(as.matrix(mat1 | mat2)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 | csr1), unname(as.matrix(mat1 | mat1)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 | emat), unname(as.matrix(mat1 | eden)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(unname(as.matrix(emat | emat)), unname(as.matrix(eden | eden)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 | csr2, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 | csc2, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 | csr1, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 | emat, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat | emat, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    
    
    expect_equal(as.matrix(csr1 & csr2), unname(as.matrix(mat1 & mat2)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 & csc2), unname(as.matrix(mat1 & mat2)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 & csr1), unname(as.matrix(mat1 & mat1)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(csr1 & emat), unname(as.matrix(mat1 & eden)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_equal(as.matrix(emat & emat), unname(as.matrix(eden & eden)))
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 & csr2, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 & csc2, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 & csr1, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(csr1 & emat, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    expect_s4_class(emat & emat, "lgRMatrix")
    expect_unmodified(csr1, csr2, csc1, csc2, emat, mat1, mat2, eden)
    

    expect_equal(as.matrix(csr1^2), as.matrix(mat1^2))
    expect_equal(as.matrix(sqrt(csr1^2)), as.matrix(sqrt(mat1^2)))
    expect_equal(as.matrix(abs(csr1)), as.matrix(abs(mat1)))
    expect_s4_class(csr1^2, "dgRMatrix")
    expect_s4_class(sqrt(csr1^2), "dgRMatrix")
    expect_s4_class(abs(csr1), "dgRMatrix")

    expect_equal(as.matrix(emat^2), as.matrix(eden^2))
    expect_equal(as.matrix(sqrt(emat^2)), as.matrix(sqrt(eden^2)))
    expect_equal(as.matrix(abs(emat)), as.matrix(abs(eden)))
    expect_s4_class(emat^2, "dgRMatrix")
    expect_s4_class(sqrt(emat^2), "dgRMatrix")
    expect_s4_class(abs(emat), "dgRMatrix")
})

test_that("Operations CSR-COO", {
    coo1 <- as.coo.matrix(csr1)
    coo2 <- as.coo.matrix(csr2)
    expect_equal(unname(as.matrix(coo1 + csr2)), unname(as.matrix(mat1 + mat2)))
    expect_equal(unname(as.matrix(coo1 + emat)), unname(as.matrix(mat1 + eden)))
    expect_equal(unname(as.matrix(emat + emat)), unname(as.matrix(eden + eden)))
    expect_s4_class(coo1 + csr2, "dgRMatrix")
    expect_s4_class(coo1 + emat, "dgRMatrix")

    expect_equal(unname(as.matrix(csr2 + coo1)), unname(as.matrix(mat2 + mat1)))
    expect_equal(unname(as.matrix(emat + coo1)), unname(as.matrix(eden + mat1)))
    expect_s4_class(csr2 + coo1, "dgRMatrix")
    expect_s4_class(emat + coo1, "dgRMatrix")

    expect_equal(unname(as.matrix(coo1 - csr2)), unname(as.matrix(mat1 - mat2)))
    expect_equal(unname(as.matrix(coo1 - emat)), unname(as.matrix(mat1 - eden)))
    expect_s4_class(coo1 - csr2, "dgRMatrix")
    expect_s4_class(coo1 - emat, "dgRMatrix")

    expect_equal(unname(as.matrix(csr2 - coo1)), unname(as.matrix(mat2 - mat1)))
    expect_equal(unname(as.matrix(emat - coo1)), unname(as.matrix(eden - mat1)))
    expect_s4_class(csr2 - coo1, "dgRMatrix")
    expect_s4_class(emat - coo1, "dgRMatrix")

    expect_equal(unname(as.matrix(coo1 * csr2)), unname(as.matrix(mat1 * mat2)))
    expect_equal(unname(as.matrix(coo1 * emat)), unname(as.matrix(mat1 * eden)))
    expect_s4_class(coo1 * csr2, "dgTMatrix")
    expect_s4_class(coo1 * emat, "dgTMatrix")

    expect_equal(unname(as.matrix(csr2 * coo1)), unname(as.matrix(mat2 * mat1)))
    expect_equal(unname(as.matrix(emat * coo1)), unname(as.matrix(eden * mat1)))
    expect_s4_class(csr2 * coo1, "dgTMatrix")
    expect_s4_class(emat * coo1, "dgTMatrix")
    
    expect_equal(unname(as.matrix(coo1 | csr2)), unname(as.matrix(mat1 | mat2)))
    expect_equal(unname(as.matrix(coo1 | emat)), unname(as.matrix(mat1 | eden)))
    expect_s4_class(coo1 | csr2, "lgRMatrix")
    expect_s4_class(coo1 | emat, "lgRMatrix")
    
    expect_equal(unname(as.matrix(csr2 | coo1)), unname(as.matrix(mat2 | mat1)))
    expect_equal(unname(as.matrix(emat | coo1)), unname(as.matrix(eden | mat1)))
    expect_s4_class(csr2 | coo1, "lgRMatrix")
    expect_s4_class(emat | coo1, "lgRMatrix")
    
    expect_equal(unname(as.matrix(coo1 & csr2)), unname(as.matrix(mat1 & mat2)))
    expect_equal(unname(as.matrix(coo1 & emat)), unname(as.matrix(mat1 & eden)))
    expect_s4_class(coo1 & csr2, "lgTMatrix")
    expect_s4_class(coo1 & emat, "lgTMatrix")
    
    expect_equal(unname(as.matrix(csr2 & coo1)), unname(as.matrix(mat2 & mat1)))
    expect_equal(unname(as.matrix(emat & coo1)), unname(as.matrix(eden & mat1)))
    expect_s4_class(csr2 & coo1, "lgTMatrix")
    expect_s4_class(emat & coo1, "lgTMatrix")


    expect_equal(unname(as.matrix(coo1^2)), unname(as.matrix(mat1^2)))
    expect_equal(unname(as.matrix(sqrt(coo1^2))), unname(as.matrix(sqrt(mat1^2))))
    expect_equal(unname(as.matrix(abs(coo1))), unname(as.matrix(abs(mat1))))
    expect_s4_class(coo1^2, "dgTMatrix")
    expect_s4_class(sqrt(coo1^2), "dgTMatrix")
    expect_s4_class(abs(coo1), "dgTMatrix")
})

test_that("Operations CSR-CSC", {
    csc1 <- as.csc.matrix(csr1)
    csc2 <- as.csc.matrix(csr2)
    expect_equal(unname(as.matrix(csc1 + csr2)), unname(as.matrix(mat1 + mat2)))
    expect_equal(unname(as.matrix(csc1 + emat)), unname(as.matrix(mat1 + eden)))
    expect_equal(unname(as.matrix(emat + emat)), unname(as.matrix(eden + eden)))
    expect_equal(unname(as.matrix(csr2 + csc1)), unname(as.matrix(mat2 + mat1)))
    expect_equal(unname(as.matrix(emat + csc1)), unname(as.matrix(eden + mat1)))

    expect_equal(unname(as.matrix(csc1 - csr2)), unname(as.matrix(mat1 - mat2)))
    expect_equal(unname(as.matrix(csc1 - emat)), unname(as.matrix(mat1 - eden)))

    expect_equal(unname(as.matrix(csr2 - csc1)), unname(as.matrix(mat2 - mat1)))
    expect_equal(unname(as.matrix(emat - csc1)), unname(as.matrix(eden - mat1)))

    expect_equal(unname(as.matrix(csc1 * csr2)), unname(as.matrix(mat1 * mat2)))
    expect_equal(unname(as.matrix(csc1 * emat)), unname(as.matrix(mat1 * eden)))

    expect_equal(unname(as.matrix(csr2 * csc1)), unname(as.matrix(mat2 * mat1)))
    expect_equal(unname(as.matrix(emat * csc1)), unname(as.matrix(eden * mat1)))
    
    expect_equal(unname(as.matrix(csc1 | csr2)), unname(as.matrix(mat1 | mat2)))
    expect_equal(unname(as.matrix(csc1 | emat)), unname(as.matrix(mat1 | eden)))
    
    expect_equal(unname(as.matrix(csr2 | csc1)), unname(as.matrix(mat2 | mat1)))
    expect_equal(unname(as.matrix(emat | csc1)), unname(as.matrix(eden | mat1)))
    
    expect_equal(unname(as.matrix(csc1 & csr2)), unname(as.matrix(mat1 & mat2)))
    expect_equal(unname(as.matrix(csc1 & emat)), unname(as.matrix(mat1 & eden)))
    
    expect_equal(unname(as.matrix(csr2 & csc1)), unname(as.matrix(mat2 & mat1)))
    expect_equal(unname(as.matrix(emat & csc1)), unname(as.matrix(eden & mat1)))
})

test_that("Types of objects", {
    set.seed(1)
    Mat <- rsparsematrix(5, 3, .4)
    Mat <- as.csr.matrix(Mat)

    expect_s4_class(Mat + as.csr.matrix(Mat, binary=TRUE), "dgRMatrix")
    expect_equal(Mat + as.csr.matrix(Mat, binary=TRUE),
                 Mat + as.csr.matrix(as.csr.matrix(Mat, binary=TRUE)))
    expect_s4_class(Mat - as.csr.matrix(Mat, binary=TRUE), "dgRMatrix")
    expect_equal(Mat - as.csr.matrix(Mat, binary=TRUE),
                 Mat - as.csr.matrix(as.csr.matrix(Mat, binary=TRUE)))
    expect_s4_class(Mat * as.csr.matrix(Mat, binary=TRUE), "dgRMatrix")
    expect_equal(Mat * as.csr.matrix(Mat, binary=TRUE),
                 Mat * as.csr.matrix(as.csr.matrix(Mat, binary=TRUE)))
    expect_s4_class(Mat | as.csr.matrix(Mat, binary=TRUE), "lgRMatrix")
    expect_equal(Mat | as.csr.matrix(Mat, binary=TRUE),
                 Mat | as.csr.matrix(as.csr.matrix(Mat, binary=TRUE)))
    expect_s4_class(Mat & as.csr.matrix(Mat, binary=TRUE), "lgRMatrix")
    expect_equal(Mat & as.csr.matrix(Mat, binary=TRUE),
                 Mat & as.csr.matrix(as.csr.matrix(Mat, binary=TRUE)))

    expect_s4_class(Mat + as.csr.matrix(Mat, logical=TRUE), "dgRMatrix")
    expect_equal(Mat + as.csr.matrix(Mat, logical=TRUE),
                 Mat + as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    expect_s4_class(Mat - as.csr.matrix(Mat, logical=TRUE), "dgRMatrix")
    expect_equal(Mat - as.csr.matrix(Mat, logical=TRUE),
                 Mat - as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    expect_s4_class(Mat * as.csr.matrix(Mat, logical=TRUE), "dgRMatrix")
    expect_equal(Mat * as.csr.matrix(Mat, logical=TRUE),
                 Mat * as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    expect_s4_class(Mat | as.csr.matrix(Mat, logical=TRUE), "lgRMatrix")
    expect_equal(Mat | as.csr.matrix(Mat, logical=TRUE),
                 Mat | as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    expect_s4_class(Mat & as.csr.matrix(Mat, logical=TRUE), "lgRMatrix")
    expect_equal(Mat & as.csr.matrix(Mat, logical=TRUE),
                 Mat & as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))




    expect_s4_class(as.csr.matrix(Mat, binary=TRUE) + as.csr.matrix(Mat, logical=TRUE), "dgRMatrix")
    expect_equal(as.csr.matrix(Mat, binary=TRUE) + as.csr.matrix(Mat, logical=TRUE),
                 as.csr.matrix(Mat, binary=TRUE) + as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    expect_s4_class(as.csr.matrix(Mat, binary=TRUE) - as.csr.matrix(Mat, logical=TRUE), "dgRMatrix")
    expect_equal(as.csr.matrix(Mat, binary=TRUE) - as.csr.matrix(Mat, logical=TRUE),
                 as.csr.matrix(Mat, binary=TRUE) - as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    expect_s4_class(as.csr.matrix(Mat, binary=TRUE) * as.csr.matrix(Mat, logical=TRUE), "dgRMatrix")
    expect_equal(as.csr.matrix(Mat, binary=TRUE) * as.csr.matrix(Mat, logical=TRUE),
                 as.csr.matrix(Mat, binary=TRUE) * as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    
    expect_s4_class(as.csr.matrix(Mat, binary=TRUE) | as.csr.matrix(Mat, logical=TRUE), "lgRMatrix")
    expect_equal(as.csr.matrix(Mat, binary=TRUE) | as.csr.matrix(Mat, logical=TRUE),
                 as.csr.matrix(Mat, binary=TRUE) | as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))
    expect_s4_class(as.csr.matrix(Mat, binary=TRUE) & as.csr.matrix(Mat, logical=TRUE), "lgRMatrix")
    expect_equal(as.csr.matrix(Mat, binary=TRUE) & as.csr.matrix(Mat, logical=TRUE),
                 as.csr.matrix(Mat, binary=TRUE) & as.csr.matrix(as.csr.matrix(Mat, logical=TRUE)))

})

test_that("CSC and dense", {
    Dense <- matrix(1:10, nrow=5L)
    Sparse <- as.csc.matrix(matrix(c(0,0,1,2,3, 11,12,13,14,15),
                                   nrow=5L, ncol=2L, byrow=FALSE))
    expect_equal(as.matrix(Dense * Sparse), Dense * as.matrix(Sparse))
    expect_s4_class(Dense * Sparse, "dgCMatrix")
    expect_equal(as.matrix(Sparse * Dense), Dense * as.matrix(Sparse))
    expect_s4_class(Sparse * Dense, "dgCMatrix")
    
    DenseNew <- Dense
    DenseNew[1,1] <- NA_real_
    
    set_new_matrix_behavior()
    expect_equal(as.matrix(DenseNew * Sparse), Dense * as.matrix(Sparse))
    expect_s4_class(DenseNew * Sparse, "dgCMatrix")
    expect_equal(as.matrix(Sparse * DenseNew), Dense * as.matrix(Sparse))
    expect_s4_class(Sparse * DenseNew, "dgCMatrix")
    
    restore_old_matrix_behavior()
    expect_equal(as.matrix(DenseNew * Sparse), DenseNew * as.matrix(Sparse))
    expect_s4_class(DenseNew * Sparse, "dgCMatrix")
    expect_equal(as.matrix(Sparse * DenseNew), DenseNew * as.matrix(Sparse))
    expect_s4_class(Sparse * DenseNew, "dgCMatrix")
    
    
    set_new_matrix_behavior()
    expect_equal(as.matrix(DenseNew & Sparse), unname(DenseNew & as.matrix(Sparse)))
    expect_s4_class(DenseNew & Sparse, "lgCMatrix")
    expect_equal(as.matrix(Sparse & DenseNew), unname(as.matrix(Sparse) & DenseNew))
    expect_s4_class(Sparse & DenseNew, "lgCMatrix")
    
    restore_old_matrix_behavior()
    expect_equal(as.matrix(DenseNew & Sparse), unname(DenseNew & as.matrix(Sparse)))
    expect_s4_class(DenseNew & Sparse, "lgCMatrix")
    expect_equal(as.matrix(Sparse & DenseNew), unname(as.matrix(Sparse) & DenseNew))
    expect_s4_class(Sparse & DenseNew, "lgCMatrix")
})

test_that("NAs in multiplication and ampersand", {
    Dense <- matrix(1:10, nrow=5L)
    Sparse <- as.csc.matrix(matrix(c(0,0,1,2,3, 11,12,13,14,15),
                                   nrow=5L, ncol=2L, byrow=FALSE))
    DenseNew <- Dense
    DenseNew[1,1] <- NA_real_
    DenseNew[2,2] <- 0
    DenseNew[3,2] <- 0
    DenseNew[5,2] <- NA_real_
    Sparse <- as.matrix(Sparse)
    Sparse[2,2] <- NA_real_
    Sparse[4,2] <- NA_real_
    Sparse[5,2] <- NA_real_
    Sparse <- as.csc.matrix(Sparse)
    DenseFilled <- DenseNew
    DenseFilled[is.na(DenseNew)] <- 0
    SparseFilled <- as.matrix(Sparse)
    SparseFilled[is.na(DenseNew) & !is.na(SparseFilled)] <- 0
    
    run_tests <- function(Sparse, DenseNew, SparseFilled, DenseFilled) {
        set_new_matrix_behavior()
        expect_equal(unname(as.matrix(DenseNew * Sparse)), unname(DenseFilled * SparseFilled))
        expect_equal(unname(as.matrix(Sparse * DenseNew)), unname(DenseFilled * SparseFilled))
        
        restore_old_matrix_behavior()
        expect_equal(unname(as.matrix(DenseNew * Sparse)), unname(DenseNew * as.matrix(Sparse)))
        expect_equal(unname(as.matrix(Sparse * DenseNew)), unname(DenseNew * as.matrix(Sparse)))
        
        set_new_matrix_behavior()
        expect_equal(as.matrix(DenseNew & Sparse), unname(DenseNew & as.matrix(Sparse)))
        expect_equal(as.matrix(Sparse & DenseNew), unname(as.matrix(Sparse) & DenseNew))
        
        restore_old_matrix_behavior()
        expect_equal(as.matrix(DenseNew & Sparse), unname(DenseNew & as.matrix(Sparse)))
        expect_equal(as.matrix(Sparse & DenseNew), unname(as.matrix(Sparse) & DenseNew))
    }
    
    run_tests(as.coo.matrix(Sparse), DenseNew, SparseFilled, DenseFilled)
    run_tests(as.csr.matrix(Sparse), DenseNew, SparseFilled, DenseFilled)
    run_tests(as.csc.matrix(Sparse), DenseNew, SparseFilled, DenseFilled)
})
