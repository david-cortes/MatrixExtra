#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <type_traits>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
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
    bool from_cpp_vec = false;
    size_t size = 0;
    std::vector<int> *int_vec_from = NULL;
    std::vector<double> *num_vec_from = NULL;
};

bool check_is_sorted(int* vec, size_t n);
bool check_is_seq(Rcpp::IntegerVector indices);
SEXP SafeRcppVector(void *args_);
bool is_same_ngRMatrix(Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
                       Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2);

/* rbind.cpp */
enum RbindedType {dgRMatrix, lgRMatrix, ngRMatrix};

/* slice.cpp */
double extract_single_val_csr
(
    int *restrict indptr,
    int *restrict indices,
    double *restrict values,
    const int row, const int col,
    const bool check_sorted
);
