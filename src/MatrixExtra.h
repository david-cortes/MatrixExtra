#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <vector>
#include <type_traits>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#ifdef _OPENMP
#   include <omp.h>
#else
#   define omp_get_thread_num() (0)
#endif

/* Aliasing for compiler optimizations */
#if defined(__GNUG__) || defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
    #define restrict __restrict
#else
    #define restrict 
#endif


#include <Rcpp.h>
#include <Rcpp/unwindProtect.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(unwindProtect)]]

extern "C" {
    #include <R.h>
    #include <Rinternals.h>
    #include <R_ext/BLAS.h>
}


/* misc.cpp */
struct VectorConstructorArgs {
    bool as_integer = false;
    bool as_logical = false;
    bool from_cpp_vec = false;
    bool from_pointer = false;
    bool cpp_lim_size = false;
    size_t size = 0;
    void *int_vec_from = NULL;
    void *num_vec_from = NULL;
    void *int_pointer_from = NULL;
    void *num_pointer_from = NULL;
};

bool check_is_sorted(int* vec, size_t n);
bool check_is_seq(Rcpp::IntegerVector indices);
SEXP SafeRcppVector(void *args_);
bool is_same_ngRMatrix(Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
                       Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2);
bool contains_any_nas_or_inf(Rcpp::NumericVector x);

/* rbind.cpp */
enum RbindedType {dgRMatrix, lgRMatrix, ngRMatrix};

/* slice.cpp */
size_t get_size_reserve(size_t nnz, size_t take1, size_t take2);
double extract_single_val_csr
(
    int *restrict indptr,
    int *restrict indices,
    double *restrict values,
    const int row, const int col,
    const bool is_sorted
);
int extract_single_val_csr
(
    int *restrict indptr,
    int *restrict indices,
    int *restrict values,
    const int row, const int col,
    const bool is_sorted
);

#if SIZE_MAX < MAX_UINT64
#   define size_large uint64_t
#else
#   define size_large size_t
#endif
