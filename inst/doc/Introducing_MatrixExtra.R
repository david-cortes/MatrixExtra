## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ---- echo=FALSE--------------------------------------------------------------
### Don't want to run compute-heavy code on CRAN
### Should build it with 'R CMD build --no-build-vignettes MatrixExtra'
### Then set this to TRUE is the vignette is ever to be rebuilt
RUN_ALL <- FALSE
if (!RUN_ALL) {
    options("MatrixExtra.nthreads" = 1L)
}

## -----------------------------------------------------------------------------
library(Matrix)
### Will construct this Matrix
### [ 1, 0, 2 ]
### [ 0, 0, 3 ]
### [ 0, 4, 0 ]
### Non-zero coordinates are:
### [(1,1), (1,3), (2,3), (3,2)]
### Row and column coordinates go separate
row_ix <- c(1, 1, 2, 3)
col_ix <- c(1, 3, 3, 2)
values <- c(1, 2, 3, 4)
X <- Matrix::sparseMatrix(
    i=row_ix, j=col_ix, x=values,
    index1=TRUE, repr="T"
)
X

## -----------------------------------------------------------------------------
as(X, "RsparseMatrix")

## -----------------------------------------------------------------------------
X + X

## -----------------------------------------------------------------------------
Xr <- as(X, "RsparseMatrix")
### This will forcibly convert the matrix to triplets
Xr[1:2, ]

## -----------------------------------------------------------------------------
### This will forcibly convert the matrix to CSC
rbind(Xr, Xr)

## -----------------------------------------------------------------------------
### This will forcibly convert the matrix to CSC
X * X

## ---- eval=RUN_ALL------------------------------------------------------------
#  library(microbenchmark)
#  set.seed(1)
#  X_big_csc <- Matrix::rsparsematrix(1e4, 1e4, .05, repr="C")
#  X_big_csr <- as(t(X_big_csc), "RsparseMatrix")
#  microbenchmark({X_slice <- X_big_csr[1:10, ]}, times=10L)

## ---- eval=RUN_ALL------------------------------------------------------------
#  microbenchmark({X_slice <- X_big_csc[, 1:10]}, times=10L)

## ---- eval=RUN_ALL------------------------------------------------------------
#  microbenchmark({X_col <- X_big_csc[, 100, drop=FALSE]}, times=10L)

## ---- eval=RUN_ALL------------------------------------------------------------
#  set.seed(1)
#  Y_dense <- matrix(rnorm(1e2*nrow(X_big_csc)), nrow=1e2)
#  microbenchmark({Z <- Y_dense %*% X_big_csc}, times=10L)

## -----------------------------------------------------------------------------
library(MatrixExtra)

## -----------------------------------------------------------------------------
Xr

## -----------------------------------------------------------------------------
options("MatrixExtra.quick_show" = FALSE)
Xr

## -----------------------------------------------------------------------------
### This will not change the format
Xr[1:2, ]

## -----------------------------------------------------------------------------
### This will not change the format
rbind(Xr, Xr)

## -----------------------------------------------------------------------------
### This will not change the format
Xr * Xr

## ---- eval=RUN_ALL------------------------------------------------------------
#  microbenchmark({X_slice <- X_big_csr[1:10, ]}, times=10L)

## ---- eval=RUN_ALL------------------------------------------------------------
#  microbenchmark({X_col <- X_big_csc[, 100, drop=FALSE]}, times=10L)

## ---- eval=RUN_ALL------------------------------------------------------------
#  microbenchmark({Z <- Y_dense %*% X_big_csc}, times=10L)

## -----------------------------------------------------------------------------
as(Xr, "ngRMatrix")

## -----------------------------------------------------------------------------
MatrixExtra::as.csr.matrix(Xr, binary=TRUE)

## -----------------------------------------------------------------------------
### Here Matrix would return a 'dgRMatrix'
t(Xr)

## -----------------------------------------------------------------------------
### Here Matrix would return a dense vector
Xr[1,]

## -----------------------------------------------------------------------------
restore_old_matrix_behavior()
set_new_matrix_behavior()

## ---- eval=RUN_ALL------------------------------------------------------------
#  library(readsparse)
#  data <- readsparse::read.sparse("real-sim")
#  X <- data$X
#  y <- as.numeric(factor(data$y))-1 ### convert to 0/1
#  X

## ---- eval=RUN_ALL------------------------------------------------------------
#  X <- cbind(rep(1, nrow(X)), X) ### Accelerated by 'MatrixExtra'
#  set.seed(1)
#  ix_train <- sample(nrow(X), floor(.5*nrow(X)), replace=FALSE)
#  X_train <- X[ix_train,] ### Accelerated by 'MatrixExtra'
#  y_train <- y[ix_train]
#  X_test <- X[-ix_train,] ### Accelerated by 'MatrixExtra'
#  y_test <- y[-ix_train]

## ---- eval=RUN_ALL------------------------------------------------------------
#  logistic_fun <- function(coefs, X, y, lambda) {
#      pred <- 1 / (1 + exp(-as.numeric(X %*% coefs))) ### Accelerated by 'MatrixExtra'
#      ll <- mean(y * log(pred) + (1 - y) * log(1 - pred))
#      reg <- lambda * as.numeric(coefs %*% coefs)
#      ### Don't regularize the intercept
#      reg <- reg - lambda * (coefs[1]^2)
#      return(-ll + reg)
#  }
#  
#  logistic_grad <- function(coefs, X, y, lambda) {
#      pred <- 1 / (1 + exp(-(X %*% coefs))) ### Accelerated by 'MatrixExtra'
#      grad <- colMeans(X * as.numeric(pred - y)) ### Accelerated by 'MatrixExtra'
#      grad <- grad + 2 * lambda * as.numeric(coefs)
#      ### Don't regularize the intercept
#      grad[1] <- grad[1] - 2 * lambda * coefs[1]
#      return(as.numeric(grad))
#  }
#  
#  lambda <- 1e-5 ### <- Regularization parameter
#  res <- optim(numeric(ncol(X_train)),
#               logistic_fun,
#               logistic_grad,
#               method="L-BFGS-B",
#               X_train, y_train, lambda)
#  fitted_coefs <- res$par

## ---- eval=RUN_ALL------------------------------------------------------------
#  y_hat_test <- as.numeric(X_test %*% fitted_coefs)
#  MLmetrics::AUC(y_hat_test, y_test)

## ---- eval=RUN_ALL------------------------------------------------------------
#  x0 <- numeric(ncol(X_train))
#  microbenchmark::microbenchmark({
#      res <- optim(x0,
#                   logistic_fun,
#                   logistic_grad,
#                   method="L-BFGS-B",
#                   X_train, y_train, lambda)
#  }, times=10L)

