library("testthat")
library("MatrixExtra")
context("RsparseMatrix subsets")

nc <- 500L
nr <- 1000L
set.seed(123)
m <- Matrix::rsparsematrix(nrow=nr, ncol=nc, density=0.1)
colnames(m) <- as.character(seq_len(nc))
rownames(m) <- as.character(seq_len(nr))
m <- as(m, "RsparseMatrix")
m_csc <- as(m, "CsparseMatrix")
m_base <- as.matrix(m)

test_that("RsparseMatrix subset cols and rows", {
    expect_equal(m, m[, ])
    expect_equal(m, m[])
    expect_equal(m, m[, , ])
    expect_equal(m[1:10, 1:100], as(m_csc[1:10, 1:100], "RsparseMatrix"))
    expect_equal(m[as.character(1:10), 1:100], as(m_csc[as.character(1:10), 1:100], "RsparseMatrix"))
    expect_equal(m["10", "20", drop=FALSE], as(m_csc["10", "20", drop=FALSE], "RsparseMatrix"))
    expect_equal(m["10", "20", drop=TRUE], m_csc["10", "20", drop=TRUE])
})

test_that("RsparseMatrix subset non sequential", {
    expect_equal(m, m[, ])
    expect_equal(m, m[])
    expect_equal(m, m[, , ])
    expect_equal(m[c(5,2,1,7,4), c(5,2,1,7,4,10,100)],
                 as(m_base[c(5,2,1,7,4), c(5,2,1,7,4,10,100)], "RsparseMatrix"))
    expect_equal(m[as.character(c(5,2,1,7,4)), as.character(c(5,2,1,7,4,10,100))],
                 as(m_base[c(5,2,1,7,4), c(5,2,1,7,4,10,100)], "RsparseMatrix"))
})

test_that("RsparseMatrix subset repeated", {
    expect_equal(m, m[, ])
    expect_equal(m, m[])
    expect_equal(m, m[, , ])
    expect_equal(m[c(2,2,2,1,1,3), c(3,3,4,4,1,1,1)],
                 as(m_base[c(2,2,2,1,1,3), c(3,3,4,4,1,1,1)], "RsparseMatrix"))
    expect_equal(m[as.character(c(2,2,2,1,1,3)), as.character(c(3,3,4,4,1,1,1))],
                 as(m_base[c(2,2,2,1,1,3), c(3,3,4,4,1,1,1)], "RsparseMatrix"))
    expect_equal(m[c(5,2,1,7,4,1,5), c(5,2,1,7,4,1,10,100,5)],
                 as(m_base[c(5,2,1,7,4,1,5), c(5,2,1,7,4,1,10,100,5)], "RsparseMatrix"))
    expect_equal(m[as.character(c(5,2,1,7,4,1,5)), as.character( c(5,2,1,7,4,1,10,100,5))],
                 as(m_base[c(5,2,1,7,4,1,5),    c(5,2,1,7,4,1,10,100,5)], "RsparseMatrix"))
})

test_that("RsparseMatrix subset empty", {
    expect_equal(m[3:10, integer()],
                 as(m_csc[3:10, integer()], "RsparseMatrix"))
    expect_equal(m[c(2,2,2,1,1,3), integer()],
                 as(m_csc[c(2,2,2,1,1,3), integer()], "RsparseMatrix"))
    expect_equal(m[, integer()],
                 as(m_csc[, integer()], "RsparseMatrix"))

    expect_equal(m[character(), ],
                 as(m_csc[integer(), ], "RsparseMatrix"))
    expect_equal(m[character(), as.character(c(3,3,4,4,1,1,1))],
                 as(m_csc[integer(), c(3,3,4,4,1,1,1)], "RsparseMatrix"))
    expect_equal(m[character(), 3:10],
                 as(m_csc[integer(), 3:10], "RsparseMatrix"))

    expect_equal(m[integer(), integer()],
                 as(m_csc[integer(), integer()], "RsparseMatrix"))
    expect_equal(m[character(), character()],
                 as(m_csc[character(), character()], "RsparseMatrix"))
})

test_that("RsparseMatrix subset cols", {
    expect_true(inherits(m[, 2L], 'numeric'))
    expect_true(inherits(m[, 2L, drop=FALSE], 'RsparseMatrix'))
    expect_true(inherits(m[, 1L:2L], 'RsparseMatrix'))
    expect_equal(rownames(m[, 2L:4L]), rownames(m))
    expect_equal(colnames(m[, 2L:4L]), as.character(2L:4L) )
    expect_equal(m[, as.character(2L:4L)], m[, 2L:4L])
    expect_error(m[, 501L])
    expect_error(m[, 500L:501L])
    expect_equal(m[, -1, drop=FALSE], as(m_csc[, -1, drop=FALSE], "RsparseMatrix"))
    expect_equal(m[, -1, drop=TRUE], as(m_csc[, -1, drop=TRUE], "RsparseMatrix"))
    expect_equal(m[, -10:-1 ], as(m_csc[, -10:-1 ], "RsparseMatrix"))
})

test_that("RsparseMatrix subset rows", {
    expect_true(inherits(m[2L, ], 'numeric'))
    expect_true(inherits(m[2L, , drop=FALSE], 'RsparseMatrix'))
    expect_true(inherits(m[1L:2L, ], 'RsparseMatrix'))
    expect_equal(colnames(m[2L:4L, ]), colnames(m))
    expect_equal(rownames(m[2L:4L, ]), as.character(2L:4L) )
    expect_equal(m[as.character(2L:4L), ], m[2L:4L, ] )
    expect_error(m[1001L, ])
    expect_error(m[900L:1001L, ])
    expect_equal(m[-1, , drop=TRUE], as(m_csc[-1, , drop=TRUE], "RsparseMatrix"))
    expect_equal(m[-1, , drop=TRUE], as(m_csc[-1, , drop=TRUE], "RsparseMatrix"))
    expect_equal(m[-10:-1, ], as(m_csc[-10:-1, ], "RsparseMatrix"))
})

test_that("RsparseMatrix subset with boolean", {
    long_vec_rows <- rep(FALSE, nrow(m))
    long_vec_cols <- rep(FALSE, ncol(m))
    long_vec_rows[1L] <- TRUE
    long_vec_rows[2L] <- TRUE
    long_vec_cols[1L] <- TRUE
    long_vec_cols[2L] <- TRUE
    expect_equal(m[long_vec_rows, ], as(m_csc[long_vec_rows, ], "RsparseMatrix"))
    expect_equal(m[, long_vec_cols], as(m_csc[, long_vec_cols], "RsparseMatrix"))
    expect_equal(m[c(TRUE, FALSE, TRUE), ], as(m_csc[c(TRUE, FALSE, TRUE), ], "RsparseMatrix"))
    expect_equal(m[, c(TRUE, FALSE, TRUE)], as(m_csc[, c(TRUE, FALSE, TRUE)], "RsparseMatrix"))
    expect_equal(m[as(c(TRUE, FALSE, TRUE), "nsparseVector"), ], as(m_csc[c(TRUE, FALSE, TRUE), ], "RsparseMatrix"))
    expect_equal(m[, as(c(TRUE, FALSE, TRUE), "nsparseVector")], as(m_csc[, c(TRUE, FALSE, TRUE)], "RsparseMatrix"))

    expect_equal(m[FALSE, ], as(m_csc[FALSE, ], "RsparseMatrix"))
    expect_equal(m[, FALSE], as(m_csc[, FALSE], "RsparseMatrix"))
    expect_equal(m[FALSE, FALSE], as(m_csc[FALSE, FALSE], "RsparseMatrix"))
    expect_equal(m[TRUE, TRUE], as(m_csc[TRUE, TRUE], "RsparseMatrix"))
})

test_that("RsparseMatrix other classes", {
    sy <- sparseMatrix(i= c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7),
                       symmetric=TRUE, dimnames = list(NULL, letters[1:7]))
    ex_dsCMatrix <- sy
    ex_lsCMatrix <- as(sy, "lsparseMatrix")
    ex_nsCMatrix <- as(sy, "nsparseMatrix")

    ex_dsRMatrix <- as(sy, "RsparseMatrix")
    ex_lsRMatrix <- as(ex_lsCMatrix, "RsparseMatrix")
    ex_nsRMatrix <- as(ex_nsCMatrix, "RsparseMatrix")

    tri <- matrix(c(1,2,0,4, 0,0,6,7, 0,0,8,9, 0,0,0,0), byrow=TRUE, nrow=4)
    tri <- as(tri, "triangularMatrix")

    ex_dtCMatrix <- as(tri, "CsparseMatrix")
    ex_ltCMatrix <- as(ex_dtCMatrix, "lsparseMatrix")
    ex_ntCMatrix <- as(ex_dtCMatrix, "nsparseMatrix")

    ex_dtRMatrix <- as(ex_dtCMatrix, "RsparseMatrix")
    ex_ltRMatrix <- as(ex_ltCMatrix, "RsparseMatrix")
    ex_ntRMatrix <- as(ex_ntCMatrix, "RsparseMatrix")

    ### Check just in case
    expect_s4_class(ex_dsRMatrix, "dsRMatrix")
    expect_s4_class(ex_lsRMatrix, "lsRMatrix")
    expect_s4_class(ex_nsRMatrix, "nsRMatrix")
    expect_s4_class(ex_dtRMatrix, "dtRMatrix")
    expect_s4_class(ex_ltRMatrix, "ltRMatrix")
    expect_s4_class(ex_ntRMatrix, "ntRMatrix")

    as.dense.matrix <- function(x) {
        if (inherits(x, "sparseMatrix"))
            x <- as.csc.matrix(x)
        x <- as.matrix(x)
        mode(x) <- "double"
        x <- unname(as.matrix(x))
        return(x)
    }

    lst_inputs <- list(
        ex_dsRMatrix, ex_lsRMatrix, ex_nsRMatrix,
        ex_dtRMatrix, ex_ltRMatrix, ex_ntRMatrix
    )
    for (inp in lst_inputs) {
        inp_dense <- as.dense.matrix(inp)

        slice_rowseq <- inp[1:3, ]
        slice_nonseq <- inp[c(2,1,3), ]
        slice_rowcol_seq <- inp[1:3, 2:4]
        slice_rowseq_randcols <- inp[1:3, c(3,2,4)]
        slice_rand <- inp[c(2,1,3), c(3,2,4)]

        dense_rowseq <- inp_dense[1:3, ]
        dense_nonseq <- inp_dense[c(2,1,3), ]
        dense_rowcol_seq <- inp_dense[1:3, 2:4]
        dense_rowseq_randcols <- inp_dense[1:3, c(3,2,4)]
        dense_rand <- inp_dense[c(2,1,3), c(3,2,4)]

        expect_s4_class(slice_rowseq, "RsparseMatrix")
        expect_s4_class(slice_nonseq, "RsparseMatrix")
        expect_s4_class(slice_rowcol_seq, "RsparseMatrix")
        expect_s4_class(slice_rowseq_randcols, "RsparseMatrix")
        expect_s4_class(slice_rand, "RsparseMatrix")

        expect_equal(as.dense.matrix(slice_rowseq), dense_rowseq)
        expect_equal(as.dense.matrix(slice_nonseq), dense_nonseq)
        expect_equal(as.dense.matrix(slice_rowcol_seq), dense_rowcol_seq)
        expect_equal(as.dense.matrix(slice_rowseq_randcols), dense_rowseq_randcols)
        expect_equal(as.dense.matrix(slice_rand), dense_rand)
    }
})
