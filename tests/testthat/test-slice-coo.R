library("testthat")
library("MatrixExtra")
restore_old_matrix_behavior()
context("TsparseMatrix subsets")

nc <- 500L
nr <- 1000L
set.seed(123)
m <- Matrix::rsparsematrix(nrow=nr, ncol=nc, density=0.1)
colnames(m) <- as.character(seq_len(nc))
rownames(m) <- as.character(seq_len(nr))
m_csc <- as(m, "CsparseMatrix")
m_coo <- as(m, "TsparseMatrix")
m_base <- as.matrix(m)
rm(m)

test_that("TsparseMatrix subset cols and rows", {
    expect_equal(m_coo, m_coo[, ])
    expect_equal(m_coo, m_coo[])
    expect_error(m_coo[, , ])
    expect_equal(as.matrix(m_coo[1:10, 1:100]),
                 m_base[1:10, 1:100])
    expect_equal(as.matrix(m_coo[as.character(1:10), 1:100]),
                 m_base[as.character(1:10), 1:100])
    expect_equal(as.matrix(m_coo["10", "20", drop=FALSE]),
                 m_base["10", "20", drop=FALSE])
    expect_equal(m_coo["10", "20", drop=TRUE],
                 m_base["10", "20", drop=TRUE])
    expect_equal(as.matrix(m_coo[10, "20", drop=FALSE]),
                 m_base[10, "20", drop=FALSE])
    expect_equal(m_coo["10", 20, drop=TRUE],
                 m_base["10", 20, drop=TRUE])
    expect_equal(as.matrix(m_coo["10", "20", drop=FALSE]),
                 m_base["10", "20", drop=FALSE])
    expect_equal(m_coo[10, 20, drop=TRUE],
                 m_base[10, 20, drop=TRUE])


    expect_equal(as.matrix(m_coo["1000", "2", drop=FALSE]),
                 m_base["1000", "2", drop=FALSE])
    expect_equal(m_coo["1000", "2", drop=TRUE],
                 m_base["1000", "2", drop=TRUE])
    expect_equal(as.matrix(m_coo[1000, "2", drop=FALSE]),
                 m_base[1000, "2", drop=FALSE])
    expect_equal(m_coo["1000", 2, drop=TRUE],
                 m_base["1000", 2, drop=TRUE])
    expect_equal(as.matrix(m_coo["1000", "2", drop=FALSE]),
                 m_base["1000", "2", drop=FALSE])
    expect_equal(m_coo[1000, 2, drop=TRUE],
                 m_base[1000, 2, drop=TRUE])

    v1 <- m_coo[,1,drop=FALSE]
    v2 <- m_coo[1,,drop=FALSE]
    expect_s4_class(v1, "dgTMatrix")
    expect_s4_class(v2, "dgTMatrix")
    expect_true(typeof(v1[,,drop=TRUE]) == "double")
    expect_true(typeof(v2[,,drop=TRUE]) == "double")
})

test_that("TsparseMatrix subset non sequential", {
    expect_equal(m_coo, m_coo[, ])
    expect_equal(m_coo, m_coo[])
    expect_error(m_coo[, , ])
    expect_equal(as.matrix(m_coo[c(5,2,1,7,4), c(5,2,1,7,4,10,100)]),
                 m_base[c(5,2,1,7,4), c(5,2,1,7,4,10,100)])
    expect_equal(as.matrix(m_coo[as.character(c(5,2,1,7,4)), as.character(c(5,2,1,7,4,10,100))]),
                 m_base[c(5,2,1,7,4), c(5,2,1,7,4,10,100)])
})

test_that("TsparseMatrix subset repeated", {
    expect_equal(as.matrix(m_coo[c(2,2,2,1,1,3), c(3,3,4,4,1,1,1)]),
                 m_base[c(2,2,2,1,1,3), c(3,3,4,4,1,1,1)])
    expect_equal(as.matrix(m_coo[as.character(c(2,2,2,1,1,3)), as.character(c(3,3,4,4,1,1,1))]),
                 m_base[c(2,2,2,1,1,3), c(3,3,4,4,1,1,1)])
    expect_equal(as.matrix(m_coo[c(5,2,1,7,4,1,5), c(5,2,1,7,4,1,10,100,5)]),
                 m_base[c(5,2,1,7,4,1,5), c(5,2,1,7,4,1,10,100,5)])
    expect_equal(as.matrix(m_coo[as.character(c(5,2,1,7,4,1,5)), as.character( c(5,2,1,7,4,1,10,100,5))]),
                 m_base[c(5,2,1,7,4,1,5),    c(5,2,1,7,4,1,10,100,5)])
})

test_that("TsparseMatrix subset empty", {
    expect_equal(as.matrix(m_coo[3:10, integer()]),
                 m_base[3:10, integer()])
    expect_equal(as.matrix(m_coo[c(2,2,2,1,1,3), integer()]),
                 m_base[c(2,2,2,1,1,3), integer()])
    expect_equal(as.matrix(m_coo[, integer()]),
                 m_base[, integer()])

    expect_equal(as.matrix(m_coo[character(), ]),
                 m_base[integer(), ])
    expect_equal(as.matrix(m_coo[character(), as.character(c(3,3,4,4,1,1,1))]),
                 m_base[integer(), c(3,3,4,4,1,1,1)])
    expect_equal(as.matrix(m_coo[character(), 3:10]),
                 m_base[integer(), 3:10])

    expect_equal(as.matrix(m_coo[integer(), integer()]),
                 unname(m_base[integer(), integer()]))
    expect_equal(as.matrix(m_coo[character(), character()]),
                 unname(m_base[character(), character()]))
})

test_that("TsparseMatrix subset cols", {
    expect_true(inherits(m_coo[, 2L], 'numeric'))
    expect_true(inherits(m_coo[, 2L, drop=FALSE], 'TsparseMatrix'))
    expect_true(inherits(m_coo[, 1L:2L], 'TsparseMatrix'))
    expect_equal(rownames(m_coo[, 2L:4L]), rownames(m_base))
    expect_equal(colnames(m_coo[, 2L:4L]), as.character(2L:4L) )
    expect_equal(m_coo[, as.character(2L:4L)], m_coo[, 2L:4L])
    expect_error(m_coo[, 501L])
    expect_error(m_coo[, 500L:501L])
    expect_equal(as.matrix(m_coo[, -1, drop=FALSE]), m_base[, -1, drop=FALSE])
    expect_equal(as.matrix(m_coo[, -1, drop=TRUE]), m_base[, -1, drop=TRUE])
    expect_equal(as.matrix(m_coo[, -10:-1 ]), m_base[, -10:-1 ])
})

test_that("TsparseMatrix subset rows", {
    expect_true(inherits(m_coo[2L, ], 'numeric'))
    expect_true(inherits(m_coo[2L, , drop=FALSE], 'TsparseMatrix'))
    expect_true(inherits(m_coo[1L:2L, ], 'TsparseMatrix'))
    expect_equal(colnames(m_coo[2L:4L, ]), colnames(m_coo))
    expect_equal(rownames(m_coo[2L:4L, ]), as.character(2L:4L) )
    expect_equal(m_coo[as.character(2L:4L), ], m_coo[2L:4L, ] )
    expect_error(m_coo[1001L, ])
    expect_error(m_coo[900L:1001L, ])
    expect_equal(as.matrix(m_coo[-1, , drop=TRUE]), m_base[-1, , drop=TRUE])
    expect_equal(as.matrix(m_coo[-1, , drop=TRUE]), m_base[-1, , drop=TRUE])
    expect_equal(as.matrix(m_coo[-10:-1, ]), m_base[-10:-1, ])
})

test_that("TsparseMatrix subset with boolean", {
    long_vec_rows <- rep(FALSE, nrow(m_coo))
    long_vec_cols <- rep(FALSE, ncol(m_coo))
    long_vec_rows[1L] <- TRUE
    long_vec_rows[2L] <- TRUE
    long_vec_cols[1L] <- TRUE
    long_vec_cols[2L] <- TRUE
    expect_equal(as.matrix(m_coo[long_vec_rows, ]),
                 m_base[long_vec_rows, ])
    expect_equal(as.matrix(m_coo[, long_vec_cols]),
                 m_base[, long_vec_cols])
    expect_equal(as.matrix(m_coo[c(TRUE, FALSE, TRUE), ]),
                 m_base[c(TRUE, FALSE, TRUE), ])
    expect_equal(as.matrix(m_coo[, c(TRUE, FALSE, TRUE)]),
                 m_base[, c(TRUE, FALSE, TRUE)])
    expect_equal(as.matrix(m_coo[as(c(TRUE, FALSE, TRUE), "nsparseVector"), ]),
                 m_base[c(TRUE, FALSE, TRUE), ])
    expect_equal(as.matrix(m_coo[, as(c(TRUE, FALSE, TRUE), "nsparseVector")]),
                 m_base[, c(TRUE, FALSE, TRUE)])

    expect_equal(as.matrix(m_coo[FALSE, ]), m_base[FALSE, ])
    expect_equal(as.matrix(m_coo[, FALSE]), m_base[, FALSE])
    expect_equal(as.matrix(m_coo[FALSE, FALSE]), unname(m_base[FALSE, FALSE]))
    expect_equal(as.matrix(m_coo[TRUE, TRUE]), m_base[TRUE, TRUE])
})

test_that("TsparseMatrix other classes", {
    sy <- sparseMatrix(i= c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7),
                       symmetric=TRUE, dimnames = list(NULL, letters[1:7]))
    ex_dsCMatrix <- sy
    ex_lsCMatrix <- as(sy, "lsparseMatrix")
    ex_nsCMatrix <- as(sy, "nsparseMatrix")

    ex_dsTMatrix <- as(sy, "TsparseMatrix")
    ex_lsTMatrix <- as(ex_lsCMatrix, "TsparseMatrix")
    ex_nsTMatrix <- as(ex_nsCMatrix, "TsparseMatrix")

    tri <- matrix(c(1,2,0,4, 0,0,6,7, 0,0,8,9, 0,0,0,0), byrow=TRUE, nrow=4)
    tri <- as(tri, "triangularMatrix")

    ex_dtCMatrix <- as(tri, "CsparseMatrix")
    ex_ltCMatrix <- as(ex_dtCMatrix, "lsparseMatrix")
    ex_ntCMatrix <- as(ex_dtCMatrix, "nsparseMatrix")

    ex_dtTMatrix <- as(ex_dtCMatrix, "TsparseMatrix")
    ex_ltTMatrix <- as(ex_ltCMatrix, "TsparseMatrix")
    ex_ntTMatrix <- as(ex_ntCMatrix, "TsparseMatrix")

    ### Check just in case
    expect_s4_class(ex_dsTMatrix, "dsTMatrix")
    expect_s4_class(ex_lsTMatrix, "lsTMatrix")
    expect_s4_class(ex_nsTMatrix, "nsTMatrix")
    expect_s4_class(ex_dtTMatrix, "dtTMatrix")
    expect_s4_class(ex_ltTMatrix, "ltTMatrix")
    expect_s4_class(ex_ntTMatrix, "ntTMatrix")

    as.dense.matrix <- function(x) {
        x_is_numeric <- inherits(x, c("dsparseMatrix", "dsparseVector"))
        x_is_logical <- inherits(x, c("lsparseMatrix", "lsparseVector",
                                      "nsparseMatrix", "nsparseVector"))
        if (inherits(x, "sparseMatrix"))
            x <- as.csc.matrix(x)
        x <- as.matrix(x)
        if (x_is_numeric)
            mode(x) <- "double"
        else
            mode(x) <- "logical"
        x <- unname(as.matrix(x))
        return(x)
    }

    lst_inputs <- list(
        ex_dsTMatrix, ex_lsTMatrix, ex_nsTMatrix,
        ex_dtTMatrix, ex_ltTMatrix, ex_ntTMatrix
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

        expect_s4_class(slice_rowseq, "TsparseMatrix")
        expect_s4_class(slice_nonseq, "TsparseMatrix")
        expect_s4_class(slice_rowcol_seq, "TsparseMatrix")
        expect_s4_class(slice_rowseq_randcols, "TsparseMatrix")
        expect_s4_class(slice_rand, "TsparseMatrix")

        expect_equal(as.dense.matrix(slice_rowseq), dense_rowseq)
        expect_equal(as.dense.matrix(slice_nonseq), dense_nonseq)
        expect_equal(as.dense.matrix(slice_rowcol_seq), dense_rowcol_seq)
        expect_equal(as.dense.matrix(slice_rowseq_randcols), dense_rowseq_randcols)
        expect_equal(as.dense.matrix(slice_rand), dense_rand)
        
        if (inherits(inp, "sparseMatrix") && nrow(inp) >= 3 && ncol(inp) >= 3)
            expect_equal(inp[3,3,drop=TRUE], inp_dense[3,3,drop=TRUE])
        if (inherits(inp, "sparseMatrix") && nrow(inp) >= 6 && ncol(inp) >= 5)
            expect_equal(inp[6,5,drop=TRUE], inp_dense[6,5,drop=TRUE])
    }
})

test_that("Reverse sequences", {
    expect_equal(as.matrix(m_coo[rev(5:100),]),
                 m_base[rev(5:100),])
    expect_equal(as.matrix(m_coo[,rev(5:100)]),
                 m_base[,rev(5:100)])
    expect_equal(as.matrix(m_coo[rev(5:100),rev(5:100)]),
                 m_base[rev(5:100),rev(5:100)])
    expect_equal(as.matrix(m_coo[rev(5:100),3,drop=FALSE]),
                 m_base[rev(5:100),3,drop=FALSE])
    expect_equal(as.matrix(m_coo[rev(5:100),c(5,3,4)]),
                 m_base[rev(5:100),c(5,3,4)])
    expect_equal(as.matrix(m_coo[c(5,3,4),rev(5:100)]),
                 m_base[c(5,3,4),rev(5:100)])
    
    expect_equal(as.matrix(m_coo[rev(1:nrow(m_coo)),]),
                 m_base[rev(1:nrow(m_base)),])
    expect_equal(as.matrix(m_coo[,rev(1:ncol(m_coo))]),
                 m_base[,rev(1:ncol(m_base))])
    expect_equal(as.matrix(m_coo[rev(1:nrow(m_coo)),rev(1:ncol(m_coo))]),
                 m_base[rev(1:nrow(m_base)),rev(1:ncol(m_base))])
    expect_equal(as.matrix(m_coo[rev(1:nrow(m_coo)),3,drop=FALSE]),
                 m_base[rev(1:nrow(m_base)),3,drop=FALSE])
    expect_equal(as.matrix(m_coo[rev(1:nrow(m_coo)),c(5,3,4)]),
                 m_base[rev(1:nrow(m_base)),c(5,3,4)])
    expect_equal(as.matrix(m_coo[c(5,3,4),rev(1:ncol(m_coo))]),
                 m_base[c(5,3,4),rev(1:ncol(m_base))])
    
    
    expect_equal(as.matrix(m_coo[rev(5:100),4:50]),
                 m_base[rev(5:100),4:50])
    expect_equal(as.matrix(m_coo[rev(1:nrow(m_coo)),4:50]),
                 m_base[rev(1:nrow(m_base)),4:50])
    
    expect_equal(as.matrix(m_coo[4:50,rev(5:100)]),
                 m_base[4:50,rev(5:100)])
    expect_equal(as.matrix(m_coo[4:50,rev(1:ncol(m_coo))]),
                 m_base[4:50,rev(1:ncol(m_base))])
})
