library("testthat")
library("Matrix")
library("MatrixExtra")
context("Methematical operators")

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

set_new_matrix_behavior()

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
        suppressWarnings({
            set_new_matrix_behavior()
            expect_equal(unname(as.matrix(DenseNew * Sparse)),
                         unname(DenseFilled * SparseFilled))
            expect_equal(unname(as.matrix(Sparse * DenseNew)),
                         unname(DenseFilled * SparseFilled))
            
            restore_old_matrix_behavior()
            expect_equal(unname(as.matrix(DenseNew * Sparse)),
                         unname(DenseNew * as.matrix(Sparse)))
            expect_equal(unname(as.matrix(Sparse * DenseNew)),
                         unname(DenseNew * as.matrix(Sparse)))
            
            set_new_matrix_behavior()
            expect_equal(as.matrix(DenseNew & Sparse),
                         unname(DenseNew & as.matrix(Sparse)))
            expect_equal(as.matrix(Sparse & DenseNew),
                         unname(as.matrix(Sparse) & DenseNew))
            
            restore_old_matrix_behavior()
            expect_equal(as.matrix(DenseNew & Sparse),
                         unname(DenseNew & as.matrix(Sparse)))
            expect_equal(as.matrix(Sparse & DenseNew),
                         unname(as.matrix(Sparse) & DenseNew))
        })
    }
    
    run_tests(as.coo.matrix(Sparse), DenseNew, SparseFilled, DenseFilled)
    run_tests(as.csr.matrix(Sparse), DenseNew, SparseFilled, DenseFilled)
    run_tests(as.csc.matrix(Sparse), DenseNew, SparseFilled, DenseFilled)
})

test_that("Vector recycling", {
    set.seed(1)
    X <- as.csr.matrix(rsparsematrix(100, 30, .1))
    v_exact <- rnorm(100)
    v_factor <- rnorm(20)
    v_factor_larger <- rnorm(500)
    v_dense <- rnorm(3000)
    v_longer <- rnorm(6000)
    v_uneven <- rnorm(80)
    v_uneven_longer <- rnorm(111)
    v_uneven_longer2 <- rnorm(333)
    v_scalar <- rnorm(1)
    X_dense <- as.matrix(X)
    
    add_rare <- function(v, pct) {
        n <- ceiling(length(v) * pct/4)
        v[sample(length(v), n, replace=FALSE)] <-  .5
        v[sample(length(v), n, replace=FALSE)] <- (-.5)
        v[sample(length(v), n, replace=FALSE)] <- (1/3)
        v[sample(length(v), n, replace=FALSE)] <- (-1/3)
        return(v)
    }
    
    add_zeros <- function(v, pct) {
        n <- ceiling(length(v) * pct)
        v[sample(length(v), n, replace=FALSE)] <- 0
        return(v)
    }
    
    add_ones <- function(v, pct) {
        n <- ceiling(length(v) * pct/2)
        v[sample(length(v), n, replace=FALSE)] <- 1
        v[sample(length(v), n, replace=FALSE)] <- (-1)
        return(v)
    }
    
    add_inf <- function(v, pct) {
        n <- ceiling(length(v) * pct/2)
        v[sample(length(v), n, replace=FALSE)] <-  Inf
        v[sample(length(v), n, replace=FALSE)] <- -Inf
        return(v)
    }
    
    as.lmatrix <- function(X) {
        X <- as.matrix(X)
        mode(X) <- "logical"
        return(X)
    }
    
    run_tests_ <- function(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
                           v_longer, v_uneven, v_uneven_longer,
                           v_uneven_longer2, v_scalar) {
        suppressWarnings({
            expect_equal(unname(as.matrix(X * v_exact)),
                         unname(X_dense * as.numeric(v_exact)))
            expect_equal(unname(as.matrix(X * v_factor)),
                         unname(X_dense * as.numeric(v_factor)))
            expect_equal(unname(as.matrix(X * v_factor_larger)),
                         unname(X_dense * as.numeric(v_factor_larger)))
            expect_equal(unname(as.matrix(X * v_dense)),
                         unname(X_dense * as.numeric(v_dense)))
            expect_error(X * v_longer)
            expect_error(X_dense * v_longer)
            expect_equal(unname(as.matrix(X * v_uneven)),
                         unname(X_dense * as.numeric(v_uneven)))
            expect_equal(unname(as.matrix(X * v_uneven_longer)),
                         unname(X_dense * as.numeric(v_uneven_longer)))
            expect_equal(unname(as.matrix(X * v_uneven_longer2)),
                         unname(X_dense * as.numeric(v_uneven_longer2)))
            expect_equal(unname(as.matrix(X * v_scalar)),
                         unname(X_dense * as.numeric(v_scalar)))
            
            
            expect_equal(unname(as.lmatrix(X & v_exact)),
                         unname(X_dense & as.numeric(v_exact)))
            expect_equal(unname(as.lmatrix(X & v_factor)),
                         unname(X_dense & as.numeric(v_factor)))
            expect_equal(unname(as.lmatrix(X & v_factor_larger)),
                         unname(X_dense & as.numeric(v_factor_larger)))
            expect_equal(unname(as.lmatrix(X & v_dense)),
                         unname(X_dense & as.numeric(v_dense)))
            expect_error(X & v_longer)
            expect_error(X_dense & v_longer)
            expect_equal(unname(as.lmatrix(X & v_uneven)),
                         unname(X_dense & as.numeric(v_uneven)))
            expect_equal(unname(as.lmatrix(X & v_uneven_longer)),
                         unname(X_dense & as.numeric(v_uneven_longer)))
            expect_equal(unname(as.lmatrix(X & v_uneven_longer2)),
                         unname(X_dense & as.numeric(v_uneven_longer2)))
            expect_equal(unname(as.lmatrix(X & v_scalar)),
                         unname(X_dense & as.numeric(v_scalar)))


            expect_equal(unname(as.matrix(X / v_exact)),
                         unname(X_dense / as.numeric(v_exact)))
            expect_equal(unname(as.matrix(X / v_factor)),
                         unname(X_dense / as.numeric(v_factor)))
            expect_equal(unname(as.matrix(X / v_factor_larger)),
                         unname(X_dense / as.numeric(v_factor_larger)))
            expect_equal(unname(as.matrix(X / v_dense)),
                         unname(X_dense / as.numeric(v_dense)))
            expect_error(X / v_longer)
            expect_error(X_dense / v_longer)
            expect_equal(unname(as.matrix(X / v_uneven)),
                         unname(X_dense / as.numeric(v_uneven)))
            expect_equal(unname(as.matrix(X / v_uneven_longer)),
                         unname(X_dense / as.numeric(v_uneven_longer)))
            expect_equal(unname(as.matrix(X / v_uneven_longer2)),
                         unname(X_dense / as.numeric(v_uneven_longer2)))
            expect_equal(unname(as.matrix(X / v_scalar)),
                         unname(X_dense / as.numeric(v_scalar)))
            
            if (!getOption("MatrixExtra.ignore_na")) {
                expect_equal(unname(as.matrix(v_exact / X)),
                             unname(as.numeric(v_exact) / X_dense))
                expect_equal(unname(as.matrix(v_factor / X)),
                             unname(as.numeric(v_factor) / X_dense))
                expect_equal(unname(as.matrix(v_factor_larger / X)),
                             unname(as.numeric(v_factor_larger) / X_dense))
                expect_equal(unname(as.matrix(v_dense / X)),
                             unname(as.numeric(v_dense) / X_dense))
                expect_error(v_longer / X)
                expect_error(v_longer / X_dense)
                expect_equal(unname(as.matrix(v_uneven / X)),
                             unname(as.numeric(v_uneven) / X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer / X)),
                             unname(as.numeric(v_uneven_longer) / X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer2 / X)),
                             unname(as.numeric(v_uneven_longer2) / X_dense))
                expect_equal(unname(as.matrix(v_scalar / X)),
                             unname(as.numeric(v_scalar) / X_dense))
            }

            if (!getOption("MatrixExtra.ignore_na")) {
                expect_equal(unname(as.matrix(X ^ v_exact)),
                             unname(X_dense ^ as.numeric(v_exact)))
                expect_equal(unname(as.matrix(X ^ v_factor)),
                             unname(X_dense ^ as.numeric(v_factor)))
                expect_equal(unname(as.matrix(X ^ v_factor_larger)),
                             unname(X_dense ^ as.numeric(v_factor_larger)))
                expect_equal(unname(as.matrix(X ^ v_dense)),
                             unname(X_dense ^ as.numeric(v_dense)))
                expect_error(X ^ v_longer)
                expect_error(X_dense ^ v_longer)
                expect_equal(unname(as.matrix(X ^ v_uneven)),
                             unname(X_dense ^ as.numeric(v_uneven)))
                expect_equal(unname(as.matrix(X ^ v_uneven_longer)),
                             unname(X_dense ^ as.numeric(v_uneven_longer)))
                expect_equal(unname(as.matrix(X ^ v_uneven_longer2)),
                             unname(X_dense ^ as.numeric(v_uneven_longer2)))
                expect_equal(unname(as.matrix(X ^ v_scalar)),
                             unname(X_dense ^ as.numeric(v_scalar)))

                expect_equal(unname(as.matrix(v_exact ^ X)),
                             unname(as.numeric(v_exact) ^ X_dense))
                expect_equal(unname(as.matrix(v_factor ^ X)),
                             unname(as.numeric(v_factor) ^ X_dense))
                expect_equal(unname(as.matrix(v_factor_larger ^ X)),
                             unname(as.numeric(v_factor_larger) ^ X_dense))
                expect_equal(unname(as.matrix(v_dense ^ X)),
                             unname(as.numeric(v_dense) ^ X_dense))
                expect_error(v_longer ^ X)
                expect_error(v_longer ^ X_dense)
                expect_equal(unname(as.matrix(v_uneven ^ X)),
                             unname(as.numeric(v_uneven) ^ X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer ^ X)),
                             unname(as.numeric(v_uneven_longer) ^ X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer2 ^ X)),
                             unname(as.numeric(v_uneven_longer2) ^ X_dense))
                expect_equal(unname(as.matrix(v_scalar ^ X)),
                             unname(as.numeric(v_scalar) ^ X_dense))
            }

            expect_equal(unname(as.matrix(X %% v_exact)),
                         unname(X_dense %% as.numeric(v_exact)))
            expect_equal(unname(as.matrix(X %% v_factor)),
                         unname(X_dense %% as.numeric(v_factor)))
            expect_equal(unname(as.matrix(X %% v_factor_larger)),
                         unname(X_dense %% as.numeric(v_factor_larger)))
            expect_equal(unname(as.matrix(X %% v_dense)),
                         unname(X_dense %% as.numeric(v_dense)))
            expect_error(X %% v_longer)
            expect_error(X_dense %% v_longer)
            expect_equal(unname(as.matrix(X %% v_uneven)),
                         unname(X_dense %% as.numeric(v_uneven)))
            expect_equal(unname(as.matrix(X %% v_uneven_longer)),
                         unname(X_dense %% as.numeric(v_uneven_longer)))
            expect_equal(unname(as.matrix(X %% v_uneven_longer2)),
                         unname(X_dense %% as.numeric(v_uneven_longer2)))
            expect_equal(unname(as.matrix(X %% v_scalar)),
                         unname(X_dense %% as.numeric(v_scalar)))
            
            if (!getOption("MatrixExtra.ignore_na")) {
                expect_equal(unname(as.matrix(v_exact %% X)),
                             unname(as.numeric(v_exact) %% X_dense))
                expect_equal(unname(as.matrix(v_factor %% X)),
                             unname(as.numeric(v_factor) %% X_dense))
                expect_equal(unname(as.matrix(v_factor_larger %% X)),
                             unname(as.numeric(v_factor_larger) %% X_dense))
                expect_equal(unname(as.matrix(v_dense %% X)),
                             unname(as.numeric(v_dense) %% X_dense))
                expect_error(v_longer %% X)
                expect_error(v_longer %% X_dense)
                expect_equal(unname(as.matrix(v_uneven %% X)),
                             unname(as.numeric(v_uneven) %% X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer %% X)),
                             unname(as.numeric(v_uneven_longer) %% X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer2 %% X)),
                             unname(as.numeric(v_uneven_longer2) %% X_dense))
                expect_equal(unname(as.matrix(v_scalar %% X)),
                             unname(as.numeric(v_scalar) %% X_dense))
            }

            expect_equal(unname(as.matrix(X %/% v_exact)),
                         unname(X_dense %/% as.numeric(v_exact)))
            expect_equal(unname(as.matrix(X %/% v_factor)),
                         unname(X_dense %/% as.numeric(v_factor)))
            expect_equal(unname(as.matrix(X %/% v_factor_larger)),
                         unname(X_dense %/% as.numeric(v_factor_larger)))
            expect_equal(unname(as.matrix(X %/% v_dense)),
                         unname(X_dense %/% as.numeric(v_dense)))
            expect_error(X %/% v_longer)
            expect_error(X_dense %/% v_longer)
            expect_equal(unname(as.matrix(X %/% v_uneven)),
                         unname(X_dense %/% as.numeric(v_uneven)))
            expect_equal(unname(as.matrix(X %/% v_uneven_longer)),
                         unname(X_dense %/% as.numeric(v_uneven_longer)))
            expect_equal(unname(as.matrix(X %/% v_uneven_longer2)),
                         unname(X_dense %/% as.numeric(v_uneven_longer2)))
            expect_equal(unname(as.matrix(X %/% v_scalar)),
                         unname(X_dense %/% as.numeric(v_scalar)))
            
            if (!getOption("MatrixExtra.ignore_na")) {
                expect_equal(unname(as.matrix(v_exact %/% X)),
                             unname(as.numeric(v_exact) %/% X_dense))
                expect_equal(unname(as.matrix(v_factor %/% X)),
                             unname(as.numeric(v_factor) %/% X_dense))
                expect_equal(unname(as.matrix(v_factor_larger %/% X)),
                             unname(as.numeric(v_factor_larger) %/% X_dense))
                expect_equal(unname(as.matrix(v_dense %/% X)),
                             unname(as.numeric(v_dense) %/% X_dense))
                expect_error(v_longer %/% X)
                expect_error(v_longer %/% X_dense)
                expect_equal(unname(as.matrix(v_uneven %/% X)),
                             unname(as.numeric(v_uneven) %/% X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer %/% X)),
                             unname(as.numeric(v_uneven_longer) %/% X_dense))
                expect_equal(unname(as.matrix(v_uneven_longer2 %/% X)),
                             unname(as.numeric(v_uneven_longer2) %/% X_dense))
                expect_equal(unname(as.matrix(v_scalar %/% X)),
                             unname(as.numeric(v_scalar) %/% X_dense))
            }
        })
    }
    
    run_tests <- function(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
                          v_longer, v_uneven, v_uneven_longer,
                          v_uneven_longer2, v_scalar) {
        run_tests_(as.csr.matrix(X), X_dense, v_exact, v_factor, v_factor_larger, v_dense,
                   v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
        run_tests_(as.coo.matrix(X), X_dense, v_exact, v_factor, v_factor_larger, v_dense,
                   v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
        
        ### Make sure that it didn't introduce any circular reference in 'Matrix'
        ### Note: 'Matrix' itself does not match exactly with base R, so these tests
        ### are bound to fail. However, they should fail due to mixing up NAs and zeros,
        ### not due to failure to execute. Thus, check them once in a while but do not
        ### leave them turned on for CRAN checks.
        # run_tests_(as.csc.matrix(X), X_dense, v_exact, v_factor, v_factor_larger, v_dense,
        #            v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
        
        # v_exact <- as(v_exact, "sparseVector")
        # v_factor <- as(v_factor, "sparseVector")
        # v_factor_larger <- as(v_factor_larger, "sparseVector")
        # v_longer <- as(v_longer, "sparseVector")
        # v_uneven <- as(v_uneven, "sparseVector")
        # v_uneven_longer <- as(v_uneven_longer, "sparseVector")
        # v_uneven_longer2 <- as(v_uneven_longer2, "sparseVector")
        # v_scalar <- as(v_scalar, "sparseVector")
        # run_tests_(as.csr.matrix(X), X_dense, v_exact, v_factor, v_factor_larger, v_dense,
        #            v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
        # run_tests_(as.csc.matrix(X), X_dense, v_exact, v_factor, v_factor_larger, v_dense,
        #            v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
        # run_tests_(as.coo.matrix(X), X_dense, v_exact, v_factor, v_factor_larger, v_dense,
        #            v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    }
    
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    v_exact <- add_rare(v_exact, .3)
    v_factor <- add_rare(v_factor, .3)
    v_factor_larger <- add_rare(v_factor_larger, .3)
    v_dense <- add_rare(v_dense, .3)
    v_uneven <- add_rare(v_uneven, .3)
    v_uneven_longer <- add_rare(v_uneven_longer, .3)
    v_uneven_longer2 <- add_rare(v_uneven_longer2, .3)
    X@x <- add_rare(X@x, .15)
    X_dense <- as.matrix(X)
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    v_exact <- add_zeros(v_exact, .3)
    v_factor <- add_zeros(v_factor, .3)
    v_factor_larger <- add_zeros(v_factor_larger, .3)
    v_dense <- add_zeros(v_dense, .3)
    v_uneven <- add_zeros(v_uneven, .3)
    v_uneven_longer <- add_zeros(v_uneven_longer, .3)
    v_uneven_longer2 <- add_zeros(v_uneven_longer2, .3)
    X@x <- add_zeros(X@x, .15)
    X_dense <- as.matrix(X)
    restore_old_matrix_behavior()
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    v_exact <- add_ones(v_exact, .15)
    v_factor <- add_ones(v_factor, .15)
    v_factor_larger <- add_ones(v_factor_larger, .15)
    v_dense <- add_ones(v_dense, .15)
    v_uneven <- add_ones(v_uneven, .15)
    v_uneven_longer <- add_ones(v_uneven_longer, .15)
    v_uneven_longer2 <- add_ones(v_uneven_longer2, .15)
    X@x <- add_ones(X@x, .1)
    X_dense <- as.matrix(X)
    restore_old_matrix_behavior()
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    v_exact[c(1,5,7,10,length(v_exact))] <- NA_real_
    v_factor[c(1,5,7,10,length(v_factor))] <- NA_real_
    v_factor_larger[c(1,5,7,10,length(v_factor_larger))] <- NA_real_
    v_dense[c(1,5,7,10,length(v_dense))] <- NA_real_
    v_longer[c(1,5,7,10,length(v_longer))] <- NA_real_
    v_uneven[c(1,5,7,10,length(v_uneven))] <- NA_real_
    v_uneven_longer[c(1,5,7,10,length(v_uneven_longer))] <- NA_real_
    v_uneven_longer2[c(1,5,7,10,length(v_uneven_longer2))] <- NA_real_
    X@x[c(1,5,7,10,length(X@x))] <- NA_real_
    X_dense <- as.matrix(X)
    restore_old_matrix_behavior()
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    v_exact <- add_inf(v_exact, .3)
    v_factor <- add_inf(v_factor, .3)
    v_factor_larger <- add_inf(v_factor_larger, .3)
    v_dense <- add_inf(v_dense, .3)
    v_uneven <- add_inf(v_uneven, .3)
    v_uneven_longer <- add_inf(v_uneven_longer, .3)
    v_uneven_longer2 <- add_inf(v_uneven_longer2, .3)
    v_scalar <- Inf
    X@x <- add_inf(X@x, .1)
    X_dense <- as.matrix(X)
    restore_old_matrix_behavior()
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    
    v_exact[5:20] <- NA_real_
    v_factor[5:20] <- NA_real_
    v_factor_larger[5:20] <- NA_real_
    v_dense[5:20] <- NA_real_
    v_longer[5:20] <- NA_real_
    v_uneven[5:20] <- NA_real_
    v_uneven_longer[5:20] <- NA_real_
    v_uneven_longer2[5:20] <- NA_real_
    v_scalar <- 1
    restore_old_matrix_behavior()
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    v_exact[] <- NA_real_
    v_factor[] <- NA_real_
    v_factor_larger[] <- NA_real_
    v_dense[] <- NA_real_
    v_longer[] <- NA_real_
    v_uneven[] <- NA_real_
    v_uneven_longer[] <- NA_real_
    v_uneven_longer2[] <- NA_real_
    v_scalar <- NA_real_
    X@x[] <- NA_real_
    X_dense <- as.matrix(X)
    restore_old_matrix_behavior()
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    v_exact[] <- 0
    v_factor[] <- 0
    v_factor_larger[] <- 0
    v_dense[] <- 0
    v_longer[] <- 0
    v_uneven[] <- 0
    v_uneven_longer[] <- 0
    v_uneven_longer2[] <- 0
    v_scalar <- 0
    X@x[] <- 0
    X_dense <- as.matrix(X)
    restore_old_matrix_behavior()
    run_tests(X, X_dense, v_exact, v_factor, v_factor_larger, v_dense,
              v_longer, v_uneven, v_uneven_longer, v_uneven_longer2, v_scalar)
    
    test_scalar <- function(X, X_dense, v_scalar) {
        suppressWarnings({
            expect_equal(unname(as.matrix(X * v_scalar)),
                         unname(X_dense * as.numeric(v_scalar)))
            expect_equal(unname(as.lmatrix(X & v_scalar)),
                         unname(X_dense & as.numeric(v_scalar)))
            expect_equal(unname(as.matrix(X / v_scalar)),
                         unname(X_dense / as.numeric(v_scalar)))
            expect_equal(unname(as.matrix(v_scalar / X)),
                         unname(as.numeric(v_scalar) / X_dense))
            expect_equal(unname(as.matrix(X ^ v_scalar)),
                         unname(X_dense ^ as.numeric(v_scalar)))
            expect_equal(unname(as.matrix(v_scalar ^ X)),
                         unname(as.numeric(v_scalar) ^ X_dense))
            expect_equal(unname(as.matrix(X %% v_scalar)),
                         unname(X_dense %% as.numeric(v_scalar)))
            expect_equal(unname(as.matrix(v_scalar %% X)),
                         unname(as.numeric(v_scalar) %% X_dense))
            expect_equal(unname(as.matrix(X %/% v_scalar)),
                         unname(X_dense %/% as.numeric(v_scalar)))
            expect_equal(unname(as.matrix(v_scalar %/% X)),
                         unname(as.numeric(v_scalar) %/% X_dense))
        })
    }
    
    for (v_scalar in c(-Inf, -3, -2, -1.5, -1, -.5, -1/3, 0,
                       1/3, .5, 1, 1.5, 2, 3, Inf, NA, numeric())) {
        test_scalar(as.csr.matrix(X), X_dense, v_scalar)
        test_scalar(as.coo.matrix(X), X_dense, v_scalar)
    }
})
