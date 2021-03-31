#include "MatrixExtra.h"

SEXP SafeRcppVector(void *args_)
{
    VectorConstructorArgs *args = (VectorConstructorArgs*)args_;
    if (args->as_integer) {
        if (args->from_cpp_vec) {
            return Rcpp::IntegerVector(args->int_vec_from->begin(), args->int_vec_from->end());
        }

        else {
            return Rcpp::IntegerVector(args->size);
        }
    }

    else {
        if (args->from_cpp_vec) {
            return Rcpp::NumericVector(args->num_vec_from->begin(), args->num_vec_from->end());
        }

        else {
            return Rcpp::NumericVector(args->size);
        }
    }
}

// [[Rcpp::export(rng = false)]]
bool is_same_ngRMatrix(Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
                       Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2)
{
  return indptr1.size() == indptr2.size() &&
         indices1.size() == indices2.size() &&
         INTEGER(indptr1) == INTEGER(indptr2) &&
         INTEGER(indices1) == INTEGER(indices2);
}

bool check_is_sorted(int* vec, size_t n)
{
    if (n <= 1)
        return true;
    if (vec[n-1] < vec[0])
        return false;
    for (size_t ix = 1; ix < n; ix++)
        if (vec[ix] < vec[ix-1])
            return false;
    return true;
}

template <class T>
void sort_sparse_indices
(
    int *restrict indptr,
    int *restrict indices,
    T values[],
    size_t nrows
)
{
    std::vector<size_t> argsorted;
    std::vector<int> temp_indices;
    std::vector<T> temp_values;
    size_t ix1, ix2;
    size_t n_this;

    for (size_t ix = 1; ix <= nrows; ix++)
    {
        ix1 = indptr[ix-1];
        ix2 = indptr[ix];
        n_this = ix2 - ix1;
        if (n_this)
        {
            if (!check_is_sorted(indices + ix1, n_this))
            {
                if (argsorted.size() < n_this) {
                    argsorted.resize(n_this);
                    temp_indices.resize(n_this);
                    temp_values.resize(n_this);
                }
                std::iota(argsorted.begin(), argsorted.begin() + n_this, (size_t)ix1);
                std::sort(argsorted.begin(), argsorted.begin() + n_this,
                          [&indices](const size_t a, const size_t b){return indices[a] < indices[b];});
                for (size_t ix = 0; ix < n_this; ix++)
                    temp_indices[ix] = indices[argsorted[ix]];
                std::copy(temp_indices.begin(), temp_indices.begin() + n_this, indices + ix1);
                for (size_t ix = 0; ix < n_this; ix++)
                    temp_values[ix] = values[argsorted[ix]];
                std::copy(temp_values.begin(), temp_values.begin() + n_this, values + ix1);

            }
        }
    }
}

void sort_sparse_indices
(
    int *restrict indptr,
    int *restrict indices,
    size_t nrows
)
{
    size_t ix1, ix2;
    size_t n_this;

    for (size_t ix = 1; ix <= nrows; ix++)
    {
        ix1 = indptr[ix-1];
        ix2 = indptr[ix];
        n_this = ix2 - ix1;
        if (n_this)
        {
            if (!check_is_sorted(indices + ix1, n_this))
            {
                std::sort(indices + ix1, indices + ix2);
            }
        }
    }
}

// [[Rcpp::export(rng = false)]]
void sort_sparse_indices_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values
)
{
    sort_sparse_indices(
        INTEGER(indptr),
        INTEGER(indices),
        REAL(values),
        indptr.size()-1
    );
}

// [[Rcpp::export(rng = false)]]
void sort_sparse_indices_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values
)
{
    sort_sparse_indices(
        INTEGER(indptr),
        INTEGER(indices),
        LOGICAL(values),
        indptr.size()-1
    );
}

// [[Rcpp::export(rng = false)]]
void sort_sparse_indices_binary
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices
)
{
    sort_sparse_indices(
        INTEGER(indptr),
        INTEGER(indices),
        indptr.size()-1
    );
}

template <class T>
void copy_from_ix(std::vector<size_t> &indices, T copy_from[], T copy_to[])
{
    for (size_t ix = 0; ix < indices.size(); ix++)
        copy_to[ix] = copy_from[indices[ix]];
}

template <class T>
void sort_coo_indices
(
    Rcpp::IntegerVector indices1,
    Rcpp::IntegerVector indices2,
    T values[]
)
{
    std::vector<size_t> argsorted(indices1.size());
    std::iota(argsorted.begin(), argsorted.end(), (size_t)0);
    int *restrict ptr1 = INTEGER(indices1);
    int *restrict ptr2 = INTEGER(indices2);
    std::sort(argsorted.begin(), argsorted.end(),
              [&ptr1, &ptr2](const size_t ix1, const size_t ix2)
              {return (ptr1[ix1] != ptr1[ix2])?
                       (ptr1[ix1] < ptr1[ix2]) : (ptr2[ix1] < ptr2[ix2]);});
    
    std::unique_ptr<char[]> temp(new char[argsorted.size() * std::max(sizeof(int), values? sizeof(T) : (size_t)1)]);
    copy_from_ix<int>(argsorted, ptr1, (int*)temp.get());
    std::copy((int*)temp.get(), (int*)temp.get() + argsorted.size(), ptr1);
    copy_from_ix<int>(argsorted, ptr2, (int*)temp.get());
    std::copy((int*)temp.get(), (int*)temp.get() + argsorted.size(), ptr2);
    if (values) {
        copy_from_ix<T>(argsorted, values, (T*)temp.get());
        std::copy((T*)temp.get(), (T*)temp.get() + argsorted.size(), values);
    }
}

// [[Rcpp::export(rng = false)]]
void sort_coo_indices_numeric
(
    Rcpp::IntegerVector indices1,
    Rcpp::IntegerVector indices2,
    Rcpp::NumericVector values
)
{
    sort_coo_indices<double>(
        indices1,
        indices2,
        REAL(values)
    );
}

// [[Rcpp::export(rng = false)]]
void sort_coo_indices_logical
(
    Rcpp::IntegerVector indices1,
    Rcpp::IntegerVector indices2,
    Rcpp::LogicalVector values
)
{
    sort_coo_indices(
        indices1,
        indices2,
        LOGICAL(values)
    );
}

// [[Rcpp::export(rng = false)]]
void sort_coo_indices_binary
(
    Rcpp::IntegerVector indices1,
    Rcpp::IntegerVector indices2
)
{
    sort_coo_indices(
        indices1,
        indices2,
        (char*)nullptr
    );
}

template <class T>
void sort_vector_indices
(
    Rcpp::IntegerVector indices_base1,
    T values[]
)
{
    std::vector<size_t> argsorted(indices_base1.size());
    std::iota(argsorted.begin(), argsorted.end(), (size_t)0);
    int *restrict ptr_indices = INTEGER(indices_base1);
    std::sort(argsorted.begin(), argsorted.end(),
              [&ptr_indices](const size_t a, const size_t b)
              {return ptr_indices[a] < ptr_indices[b];});
    std::unique_ptr<char[]> temp(new char[argsorted.size() * std::max(sizeof(int), values? sizeof(T) : (size_t)1)]);
    copy_from_ix<int>(argsorted, ptr_indices, (int*)temp.get());
    std::copy((int*)temp.get(), (int*)temp.get() + argsorted.size(), ptr_indices);
    if (values) {
        copy_from_ix<T>(argsorted, values, (T*)temp.get());
        std::copy((T*)temp.get(), (T*)temp.get() + argsorted.size(), values);
    }
}

// [[Rcpp::export(rng = false)]]
void sort_vector_indices_numeric
(
    Rcpp::IntegerVector indices_base1,
    Rcpp::NumericVector values
)
{
    sort_vector_indices(
        indices_base1,
        REAL(values)
    );
}

// [[Rcpp::export(rng = false)]]
void sort_vector_indices_integer
(
    Rcpp::IntegerVector indices_base1,
    Rcpp::IntegerVector values
)
{
    sort_vector_indices(
        indices_base1,
        INTEGER(values)
    );
}

// [[Rcpp::export(rng = false)]]
void sort_vector_indices_logical
(
    Rcpp::IntegerVector indices_base1,
    Rcpp::LogicalVector values
)
{
    sort_vector_indices(
        indices_base1,
        LOGICAL(values)
    );
}

// [[Rcpp::export(rng = false)]]
void sort_vector_indices_binary
(
    Rcpp::IntegerVector indices_base1
)
{
    std::sort(INTEGER(indices_base1), INTEGER(indices_base1) + indices_base1.size());
}
