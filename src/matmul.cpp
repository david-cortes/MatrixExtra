#include "MatrixExtra.h"

/* Mental map to figure out what should be called where:

matmul(x,y) -> x %*% y
    = ColMajor(Xc * Yc)
    = RowMajor(Yr * Xr)

crossprod(x,y) -> t(x) %*% y
    = ColMajor(Xc.T * Yc)
    = ColMajor(Xr * Yc)
    = RowMajor(Yr * Xc)

tcrossprod(x,y) -> x %*% t(y)
    = ColMajor(Xc * Yc.T)
    = ColMajor(Xc * Yr)
    = RowMajor(Yc * Xr)

The required functions are:
    matmul(CSR, dense)
    tcrossprod(CSR, dense)
    tcrossprod(dense, CSR)
    matmul(dense, CSC)
    crossprod(dense, CSC)

    with all the dense ones being column-major
*/

/* Unoptimized replacement for BLAS axpy in single precision
   
   Note: R does not ship with single-precision BLAS routines.
   Using an optimized 'saxpy' would imply having to link against
   the 'float' package, and in such case there might be problems if one
   changes the BLAS shared object within a session.
   'axpy' is anyway an operation without much room for optimizations
   beyond what the compiler already does, so the performance hit
   might be very small or non-existent. */
static inline void saxpy_cust(const size_t n, const float a,
                              const float *restrict x, const size_t incx,
                              float *restrict y, const size_t incy)
{
    #pragma omp simd
    for (size_t ix = 0; ix < n; ix++) y[ix*incy] += a * x[ix*incx];
}

/* X <- A*B + X | A(m,k) is sparse CSR, B(k,n) is dense row-major, X(m,n) is dense row-major
   
   Equivalences:
    X <- A*t(B) + X    | A(m,k) CSR, B(n,k) column-major, X(m,n) row-major
    X <- t(A)*B + X    | A(k,m) CSC, B(k,n) row-major, X(m,n) row-major
    X <- t(A)*t(B) + X | A(k,m) CSC, B(n,k) column-major, X(m,n) row-major
    
    When X is zeroed-out:
    X <- B*t(A) | A(k,m) CSR, B(n,k) column-major, X(n,m) column-major
    X <- B*A    | A(k,m) CSC, B(n,k) column-major, X(n,m) column-major

*/
void dgemm_csr_drm_as_drm
(
    const int m, const int n,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    const double *restrict DenseMat, const size_t ldb,
    double *restrict OutputMat, const size_t ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    const int one = 1;
    double *restrict row_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, ldb, ldc, OutputMat, DenseMat, indptr, indices, values, one) \
            private(row_ptr)
    for (int row = 0; row < m; row++)
    {
        row_ptr = OutputMat + (size_t)row*ldc;
        for (int col = indptr[row]; col < indptr[row+1]; col++)
            daxpy_(&n, values + col, DenseMat + (size_t)indices[col]*ldb, &one, row_ptr, &one);
    }
}

void dgemm_csr_drm_as_drm
(
    const int m, const int n,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    const float *restrict DenseMat, const size_t ldb,
    float *restrict OutputMat, const size_t ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    const int one = 1;
    float *restrict row_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, ldb, ldc, OutputMat, DenseMat, indptr, indices, values, one) \
            private(row_ptr)
    for (int row = 0; row < m; row++)
    {
        row_ptr = OutputMat + (size_t)row*ldc;
        for (int col = indptr[row]; col < indptr[row+1]; col++)
            saxpy_cust(n, values[col], DenseMat + (size_t)indices[col]*ldb, one, row_ptr, one);
    }
}

/* X <- A*B + X | A(m,k) is sparse CSR, B(k,n) is dense row-major, X(m,n) is dense column-major

   When X is zeroed-out:
   X <- A*t(B) | A(m,k) CSR, B(n,k) column-major, X(m,n) column-major

*/
void dgemm_csr_drm_as_dcm
(
    const int m, const int n,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    const double *restrict DenseMat, const size_t ldb,
    double *restrict OutputMat, const int ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    const int one = 1;
    double *restrict row_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, ldb, ldc, OutputMat, DenseMat, indptr, indices, values, one) \
            private(row_ptr)
    for (int row = 0; row < m; row++)
    {
        row_ptr = OutputMat + row;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            daxpy_(&n, values + ix, DenseMat + (size_t)indices[ix]*ldb, &one, row_ptr, &ldc);
    }
}

void dgemm_csr_drm_as_dcm
(
    const int m, const int n,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    const float *restrict DenseMat, const size_t ldb,
    float *restrict OutputMat, const int ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    const int one = 1;
    float *restrict row_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, ldb, ldc, OutputMat, DenseMat, indptr, indices, values, one) \
            private(row_ptr)
    for (int row = 0; row < m; row++)
    {
        row_ptr = OutputMat + row;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            saxpy_cust(n, values[ix], DenseMat + (size_t)indices[ix]*ldb, one, row_ptr, ldc);
    }
}

/* X <- A*B + X | A(m,k) is dense row-major, B(k,n) is sparse CSC, X(m,n) is dense column-major

   When X is zeroed-out:
   X <- t(A)*B | A(k,m) column-major, B(k,n) CSC, X(m,n) column-major

*/
void dgemm_drm_csc_as_dcm
(
    const int m, const int n,
    const double *restrict DenseMat, const int lda,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    double *restrict OutputMat, const size_t ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    const int one = 1;
    double *restrict col_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, DenseMat, lda, indptr, indices, values, OutputMat) \
            private(col_ptr)
    for (int col = 0; col < n; col++)
    {
        col_ptr = OutputMat + (size_t)col*ldc;
        for (int ix = indptr[col]; ix < indptr[col+1]; ix++)
            daxpy_(&m, values + ix, DenseMat + indices[ix], &lda, col_ptr, &one);
    }
}

void dgemm_drm_csc_as_dcm
(
    const int m, const int n,
    const float *restrict DenseMat, const int lda,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    float *restrict OutputMat, const size_t ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    const int one = 1;
    float *restrict col_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, DenseMat, lda, indptr, indices, values, OutputMat) \
            private(col_ptr)
    for (int col = 0; col < n; col++)
    {
        col_ptr = OutputMat + (size_t)col*ldc;
        for (int ix = indptr[col]; ix < indptr[col+1]; ix++)
            saxpy_cust(m, values[ix], DenseMat + indices[ix], lda, col_ptr, one);
    }
}

/* X <- A*B | A(m,k) is sparse CSR, B(k,n) is dense column-major, X(m,n) is dense column major */
void dgemm_csr_dcm_as_dcm
(
    const int m, const int n,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    const double *restrict DenseMat, const int ldb,
    double *restrict OutputMat, const int ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    double *restrict row_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, indptr, indices, values, DenseMat, ldb, OutputMat, ldc) \
            private(row_ptr)
    for (int row = 0; row < m; row++)
    {
        row_ptr = OutputMat + row;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            daxpy_(&n, values + ix, DenseMat + indices[ix], &ldb, row_ptr, &ldc);
    }
}

void dgemm_csr_dcm_as_dcm
(
    const int m, const int n,
    const int *restrict indptr, const int *restrict indices, const double *restrict values,
    const float *restrict DenseMat, const int ldb,
    float *restrict OutputMat, const int ldc,
    int nthreads
)
{
    if (m <= 0 || indptr[0] == indptr[m])
        return;
    float *restrict row_ptr;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(m, n, indptr, indices, values, DenseMat, ldb, OutputMat, ldc) \
            private(row_ptr)
    for (int row = 0; row < m; row++)
    {
        row_ptr = OutputMat + row;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            saxpy_cust(n, values[ix], DenseMat + indices[ix], ldb, row_ptr, ldc);
    }
}

/* x %*% y */
template <class RcppMatrix>
RcppMatrix matmul_dense_csc(RcppMatrix X_colmajor,
                            Rcpp::IntegerVector Y_csc_indptr,
                            Rcpp::IntegerVector Y_csc_indices,
                            Rcpp::NumericVector Y_csc_values,
                            int nthreads)
{
    int nrows_X = X_colmajor.nrow();
    int ncols_Y = Y_csc_indptr.size() - 1;
    RcppMatrix out_colmajor(nrows_X, ncols_Y);
    int nrows = out_colmajor.nrow();
    int ncols = out_colmajor.ncol();

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
        dgemm_csr_drm_as_drm(
            ncols, nrows,
            INTEGER(Y_csc_indptr), INTEGER(Y_csc_indices), REAL(Y_csc_values),
            REAL(X_colmajor), nrows,
            REAL(out_colmajor), nrows,
            nthreads
        );
    else
        dgemm_csr_drm_as_drm(
            ncols, nrows,
            INTEGER(Y_csc_indptr), INTEGER(Y_csc_indices), REAL(Y_csc_values),
            (float*)INTEGER(X_colmajor), nrows,
            (float*)INTEGER(out_colmajor), nrows,
            nthreads
        );

    return out_colmajor;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix matmul_dense_csc_numeric(Rcpp::NumericMatrix X_colmajor,
                                             Rcpp::IntegerVector Y_csc_indptr,
                                             Rcpp::IntegerVector Y_csc_indices,
                                             Rcpp::NumericVector Y_csc_values,
                                             int nthreads)
{
    return matmul_dense_csc<Rcpp::NumericMatrix>(
        X_colmajor,
        Y_csc_indptr,
        Y_csc_indices,
        Y_csc_values,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix matmul_dense_csc_float32(Rcpp::IntegerMatrix X_colmajor,
                                             Rcpp::IntegerVector Y_csc_indptr,
                                             Rcpp::IntegerVector Y_csc_indices,
                                             Rcpp::NumericVector Y_csc_values,
                                             int nthreads)
{
    return matmul_dense_csc<Rcpp::IntegerMatrix>(
        X_colmajor,
        Y_csc_indptr,
        Y_csc_indices,
        Y_csc_values,
        nthreads
    );
}

/* x %*% t(y) */
template <class RcppMatrix>
RcppMatrix tcrossprod_dense_csr(RcppMatrix X_colmajor,
                                Rcpp::IntegerVector Y_csr_indptr,
                                Rcpp::IntegerVector Y_csr_indices,
                                Rcpp::NumericVector Y_csr_values,
                                int nthreads, int ncols_Y)
{
    RcppMatrix out_colmajor(X_colmajor.nrow(), Y_csr_indptr.size()-1);

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
        dgemm_csr_drm_as_drm(
            out_colmajor.ncol(), out_colmajor.nrow(),
            INTEGER(Y_csr_indptr), INTEGER(Y_csr_indices), REAL(Y_csr_values),
            REAL(X_colmajor), X_colmajor.nrow(),
            REAL(out_colmajor), out_colmajor.nrow(),
            nthreads
        );
    else
        dgemm_csr_drm_as_drm(
            out_colmajor.ncol(), out_colmajor.nrow(),
            INTEGER(Y_csr_indptr), INTEGER(Y_csr_indices), REAL(Y_csr_values),
            (float*)INTEGER(X_colmajor), X_colmajor.nrow(),
            (float*)INTEGER(out_colmajor), out_colmajor.nrow(),
            nthreads
        );

    return out_colmajor;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix tcrossprod_dense_csr_numeric(Rcpp::NumericMatrix X_colmajor,
                                                 Rcpp::IntegerVector Y_csr_indptr,
                                                 Rcpp::IntegerVector Y_csr_indices,
                                                 Rcpp::NumericVector Y_csr_values,
                                                 int nthreads, int ncols_Y)
{
    return tcrossprod_dense_csr<Rcpp::NumericMatrix>(
        X_colmajor,
        Y_csr_indptr,
        Y_csr_indices,
        Y_csr_values,
        nthreads, ncols_Y
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix tcrossprod_dense_csr_float32(Rcpp::IntegerMatrix X_colmajor,
                                                 Rcpp::IntegerVector Y_csr_indptr,
                                                 Rcpp::IntegerVector Y_csr_indices,
                                                 Rcpp::NumericVector Y_csr_values,
                                                 int nthreads, int ncols_Y)
{
    return tcrossprod_dense_csr<Rcpp::IntegerMatrix>(
        X_colmajor,
        Y_csr_indptr,
        Y_csr_indices,
        Y_csr_values,
        nthreads, ncols_Y
    );
}

/* t(x) %*% y */
template <class RcppMatrix>
RcppMatrix crossprod_dense_csc(RcppMatrix X_colmajor,
                               Rcpp::IntegerVector Y_csc_indptr,
                               Rcpp::IntegerVector Y_csc_indices,
                               Rcpp::NumericVector Y_csc_values,
                               int nthreads)
{
    RcppMatrix out_colmajor(X_colmajor.ncol(), Y_csc_indptr.size()-1);

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
        dgemm_drm_csc_as_dcm(
            out_colmajor.nrow(), out_colmajor.ncol(),
            REAL(X_colmajor), X_colmajor.nrow(),
            INTEGER(Y_csc_indptr), INTEGER(Y_csc_indices), REAL(Y_csc_values),
            REAL(out_colmajor), out_colmajor.nrow(),
            nthreads
        );
    else
        dgemm_drm_csc_as_dcm(
            out_colmajor.nrow(), out_colmajor.ncol(),
            (float*)INTEGER(X_colmajor), X_colmajor.nrow(),
            INTEGER(Y_csc_indptr), INTEGER(Y_csc_indices), REAL(Y_csc_values),
            (float*)INTEGER(out_colmajor), out_colmajor.nrow(),
            nthreads
        );

    return out_colmajor;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix crossprod_dense_csc_numeric(Rcpp::NumericMatrix X_colmajor,
                                                Rcpp::IntegerVector Y_csc_indptr,
                                                Rcpp::IntegerVector Y_csc_indices,
                                                Rcpp::NumericVector Y_csc_values,
                                                int nthreads)
{
    return crossprod_dense_csc<Rcpp::NumericMatrix>(
        X_colmajor,
        Y_csc_indptr,
        Y_csc_indices,
        Y_csc_values,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix crossprod_dense_csc_float32(Rcpp::IntegerMatrix X_colmajor,
                                                Rcpp::IntegerVector Y_csc_indptr,
                                                Rcpp::IntegerVector Y_csc_indices,
                                                Rcpp::NumericVector Y_csc_values,
                                                int nthreads)
{
    return crossprod_dense_csc<Rcpp::IntegerMatrix>(
        X_colmajor,
        Y_csc_indptr,
        Y_csc_indices,
        Y_csc_values,
        nthreads
    );
}

/* x %*% y */
template <class RcppMatrix>
RcppMatrix matmul_csr_dense(Rcpp::IntegerVector X_csr_indptr,
                            Rcpp::IntegerVector X_csr_indices,
                            Rcpp::NumericVector X_csr_values,
                            RcppMatrix Y_colmajor,
                            int nthreads)
{
    RcppMatrix out_colmajor(X_csr_indptr.size()-1, Y_colmajor.ncol());

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
        dgemm_csr_dcm_as_dcm(
            out_colmajor.nrow(), out_colmajor.ncol(),
            INTEGER(X_csr_indptr), INTEGER(X_csr_indices), REAL(X_csr_values),
            REAL(Y_colmajor), Y_colmajor.nrow(),
            REAL(out_colmajor), out_colmajor.nrow(),
            nthreads
        );
    else
        dgemm_csr_dcm_as_dcm(
            out_colmajor.nrow(), out_colmajor.ncol(),
            INTEGER(X_csr_indptr), INTEGER(X_csr_indices), REAL(X_csr_values),
            (float*)INTEGER(Y_colmajor), Y_colmajor.nrow(),
            (float*)INTEGER(out_colmajor), out_colmajor.nrow(),
            nthreads
        );
    return out_colmajor;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix matmul_csr_dense_numeric(Rcpp::IntegerVector X_csr_indptr,
                                             Rcpp::IntegerVector X_csr_indices,
                                             Rcpp::NumericVector X_csr_values,
                                             Rcpp::NumericMatrix Y_colmajor,
                                             int nthreads)
{
    return matmul_csr_dense<Rcpp::NumericMatrix>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        Y_colmajor,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix matmul_csr_dense_float32(Rcpp::IntegerVector X_csr_indptr,
                                             Rcpp::IntegerVector X_csr_indices,
                                             Rcpp::NumericVector X_csr_values,
                                             Rcpp::IntegerMatrix Y_colmajor,
                                             int nthreads)
{
    return matmul_csr_dense<Rcpp::IntegerMatrix>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        Y_colmajor,
        nthreads
    );
}

/* x %*% t(y) */
template <class RcppMatrix>
RcppMatrix tcrossprod_csr_dense(Rcpp::IntegerVector X_csr_indptr,
                                Rcpp::IntegerVector X_csr_indices,
                                Rcpp::NumericVector X_csr_values,
                                RcppMatrix Y_colmajor,
                                int nthreads)
{
    RcppMatrix out_colmajor(X_csr_indptr.size()-1, Y_colmajor.nrow());

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
        dgemm_csr_drm_as_dcm(
            out_colmajor.nrow(), out_colmajor.ncol(),
            INTEGER(X_csr_indptr), INTEGER(X_csr_indices), REAL(X_csr_values),
            REAL(Y_colmajor), Y_colmajor.nrow(),
            REAL(out_colmajor), out_colmajor.nrow(),
            nthreads
        );
    else
        dgemm_csr_drm_as_dcm(
            out_colmajor.nrow(), out_colmajor.ncol(),
            INTEGER(X_csr_indptr), INTEGER(X_csr_indices), REAL(X_csr_values),
            (float*)INTEGER(Y_colmajor), Y_colmajor.nrow(),
            (float*)INTEGER(out_colmajor), out_colmajor.nrow(),
            nthreads
        );

    return out_colmajor;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix tcrossprod_csr_dense_numeric(Rcpp::IntegerVector X_csr_indptr,
                                                 Rcpp::IntegerVector X_csr_indices,
                                                 Rcpp::NumericVector X_csr_values,
                                                 Rcpp::NumericMatrix Y_colmajor,
                                                 int nthreads)
{
    return tcrossprod_csr_dense<Rcpp::NumericMatrix>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        Y_colmajor,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix tcrossprod_csr_dense_float32(Rcpp::IntegerVector X_csr_indptr,
                                                 Rcpp::IntegerVector X_csr_indices,
                                                 Rcpp::NumericVector X_csr_values,
                                                 Rcpp::IntegerMatrix Y_colmajor,
                                                 int nthreads)
{
    return tcrossprod_csr_dense<Rcpp::IntegerMatrix>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        Y_colmajor,
        nthreads
    );
}

/* x %*% y */
template <class RcppVector>
Rcpp::NumericVector matmul_csr_dvec(Rcpp::IntegerVector X_csr_indptr,
                                    Rcpp::IntegerVector X_csr_indices,
                                    Rcpp::NumericVector X_csr_values,
                                    RcppVector y_dense,
                                    int nthreads)
{
    Rcpp::NumericVector out(X_csr_indptr.size()-1);
    const int nrows = out.size();

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(X_csr_indptr, X_csr_indices, X_csr_values, y_dense, nrows)
    for (int row = 0; row < nrows; row++)
        for (int ix = X_csr_indptr[row]; ix < X_csr_indptr[row+1]; ix++)
        {
            if (std::is_same<RcppVector, Rcpp::IntegerVector>::value)
                out[row] += (y_dense[X_csr_indices[ix]] == NA_INTEGER)?
                             NA_REAL : X_csr_values[ix] * y_dense[X_csr_indices[ix]];
            else if (std::is_same<RcppVector, Rcpp::LogicalVector>::value)
                out[row] += (y_dense[X_csr_indices[ix]] == NA_LOGICAL)?
                             NA_REAL : X_csr_values[ix] * (bool)y_dense[X_csr_indices[ix]];
            else
                out[row] += X_csr_values[ix] * y_dense[X_csr_indices[ix]];
        }

    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_dvec_numeric(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::NumericVector y_dense,
                                            int nthreads)
{
    return matmul_csr_dvec<Rcpp::NumericVector>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_dense,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_dvec_integer(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::IntegerVector y_dense,
                                            int nthreads)
{
    return matmul_csr_dvec<Rcpp::IntegerVector>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_dense,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_dvec_logical(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::LogicalVector y_dense,
                                            int nthreads)
{
    return matmul_csr_dvec<Rcpp::LogicalVector>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_dense,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_dvec_float32(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::IntegerVector y_dense,
                                            int nthreads)
{
    return matmul_csr_dvec<float*>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        (float*)INTEGER(y_dense),
        nthreads
    );
}

/* x %*% y */
template <class RcppVector>
Rcpp::NumericVector matmul_csr_svec(Rcpp::IntegerVector X_csr_indptr,
                                    Rcpp::IntegerVector X_csr_indices,
                                    Rcpp::NumericVector X_csr_values,
                                    Rcpp::IntegerVector y_indices_base1,
                                    RcppVector y_values,
                                    int nthreads)
{
    Rcpp::NumericVector out(X_csr_indptr.size()-1);
    if (!y_indices_base1.size())
        return out;
    const int nrows = out.size();

    int *restrict ptr_X_indices = INTEGER(X_csr_indices);
    int *restrict ptr_y_indices = INTEGER(y_indices_base1);
    int *restrict end_y = ptr_y_indices + y_indices_base1.size();
    int *ptr1, *ptr2, *end1;

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            shared(nrows, X_csr_indptr, X_csr_indices, X_csr_values, \
                   ptr_X_indices, ptr_y_indices, end_y, out) \
            private(ptr1, ptr2, end1)
    for (int row = 0; row < nrows; row++)
    {
        ptr1 = ptr_X_indices + X_csr_indptr[row];
        end1 = ptr_X_indices + X_csr_indptr[row+1];
        ptr2 = ptr_y_indices;

        while (true)
        {
            if (ptr1 >= end1 || ptr2 >= end_y) {
                goto next_row;
            }

            else if (*ptr1 == (*ptr2)-1) {
                if (std::is_same<RcppVector, Rcpp::IntegerVector>::value)
                    out[row] += (y_values[ptr2 - ptr_y_indices] == NA_INTEGER)?
                                 NA_REAL : X_csr_values[ptr1 - ptr_X_indices] * y_values[ptr2 - ptr_y_indices];
                else if (std::is_same<RcppVector, Rcpp::LogicalVector>::value)
                    out[row] += (y_values[ptr2 - ptr_y_indices] == NA_LOGICAL)?
                                 NA_REAL : X_csr_values[ptr1 - ptr_X_indices] * (bool)y_values[ptr2 - ptr_y_indices];
                else if (std::is_same<RcppVector, char*>::value)
                    out[row] += X_csr_values[ptr1 - ptr_X_indices];
                else
                    out[row] += X_csr_values[ptr1 - ptr_X_indices] * y_values[ptr2 - ptr_y_indices];
                ptr1++;
                ptr2++;
            }

            else if (*ptr2-1 > *ptr1) {
                ptr1 = std::lower_bound(ptr1, end1, *ptr2-1);
            }

            else {
                ptr2 = std::lower_bound(ptr2, end_y, *ptr1+1);
            }
        }

        next_row:
        {}
    }

    return out;
}


// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_svec_numeric(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::IntegerVector y_indices_base1,
                                            Rcpp::NumericVector y_values,
                                            int nthreads)
{
    return matmul_csr_svec<Rcpp::NumericVector>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_indices_base1,
        y_values,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_svec_integer(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::IntegerVector y_indices_base1,
                                            Rcpp::IntegerVector y_values,
                                            int nthreads)
{
    return matmul_csr_svec<Rcpp::IntegerVector>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_indices_base1,
        y_values,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_svec_logical(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::IntegerVector y_indices_base1,
                                            Rcpp::LogicalVector y_values,
                                            int nthreads)
{
    return matmul_csr_svec<Rcpp::LogicalVector>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_indices_base1,
        y_values,
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_svec_binary(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::IntegerVector y_indices_base1,
                                            int nthreads)
{
    return matmul_csr_svec<char*>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_indices_base1,
        (char*)INTEGER(y_indices_base1), /* <- placeholder */
        nthreads
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector matmul_csr_svec_float32(Rcpp::IntegerVector X_csr_indptr,
                                            Rcpp::IntegerVector X_csr_indices,
                                            Rcpp::NumericVector X_csr_values,
                                            Rcpp::IntegerVector y_indices_base1,
                                            Rcpp::IntegerVector y_values,
                                            int nthreads)
{
    return matmul_csr_svec<float*>(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        y_indices_base1,
        (float*)INTEGER(y_values),
        nthreads
    );
}
