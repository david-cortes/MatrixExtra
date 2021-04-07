#include "MatrixExtra.h"

SEXP SafeRcppVector(void *args_)
{
    VectorConstructorArgs *args = (VectorConstructorArgs*)args_;
    std::vector<int> *int_vec_from = (std::vector<int>*)args->int_vec_from;
    std::vector<double> *num_vec_from = (std::vector<double>*)args->num_vec_from;
    
    if (args->as_integer) {
        if (args->from_cpp_vec) {
            if (!args->as_logical)
                return Rcpp::IntegerVector(int_vec_from->begin(), int_vec_from->end());
            else
                return Rcpp::LogicalVector(int_vec_from->begin(), int_vec_from->end());
        }

        else {
            if (!args->as_logical)
                return Rcpp::IntegerVector(args->size);
            else
                return Rcpp::LogicalVector(args->size);
        }
    }

    else {
        if (args->from_cpp_vec) {
            return Rcpp::NumericVector(num_vec_from->begin(), num_vec_from->end());
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

bool check_indices_are_unsorted
(
    int *restrict indptr,
    int *restrict indices,
    int nrows
)
{
    for (int row = 0; row < nrows; row++)
    {
        if(!check_is_sorted(indices + indptr[row], indptr[row+1] - indptr[row]))
            return false;
    }
    return true;
}

// [[Rcpp::export(rng = false)]]
bool check_indices_are_unsorted
(
    Rcpp::IntegerVector indptr,
    Rcpp::NumericVector indices
)
{
    return check_indices_are_unsorted(
        INTEGER(indptr),
        INTEGER(indices),
        indptr.size() - 1
    );
}


template <class T>
void sort_sparse_indices
(
    int *restrict indptr,
    int *restrict indices,
    T values[],
    int nrows
)
{
    std::vector<int> argsorted;
    std::vector<int> temp_indices;
    std::vector<T> temp_values;
    int ix1, ix2;
    int n_this;

    for (int row = 1; row <= nrows; row++)
    {
        ix1 = indptr[row-1];
        ix2 = indptr[row];
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
                std::iota(argsorted.begin(), argsorted.begin() + n_this, ix1);
                std::sort(argsorted.begin(), argsorted.begin() + n_this,
                          [&indices](const int a, const int b){return indices[a] < indices[b];});
                for (int ix = 0; ix < n_this; ix++)
                    temp_indices[ix] = indices[argsorted[ix]];
                std::copy(temp_indices.begin(), temp_indices.begin() + n_this, indices + ix1);
                for (int ix = 0; ix < n_this; ix++)
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

    for (size_t row = 1; row <= nrows; row++)
    {
        ix1 = indptr[row-1];
        ix2 = indptr[row];
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

template <class T>
void sort_sparse_indices_known_ncol
(
    int *restrict indptr,
    int *restrict indices,
    T values[],
    int nrows, int ncols
)
{
    std::vector<int> argsorted(ncols);
    std::vector<int> temp_indices(ncols);
    std::vector<T> temp_values(ncols);
    int ix1, ix2;
    int n_this;

    for (int row = 1; row <= nrows; row++)
    {
        ix1 = indptr[row-1];
        ix2 = indptr[row];
        n_this = ix2 - ix1;
        if (n_this)
        {
            if (!check_is_sorted(indices + ix1, n_this))
            {
                std::iota(argsorted.begin(), argsorted.begin() + n_this, ix1);
                std::sort(argsorted.begin(), argsorted.begin() + n_this,
                          [&indices](const int a, const int b){return indices[a] < indices[b];});
                for (int ix = 0; ix < n_this; ix++)
                    temp_indices[ix] = indices[argsorted[ix]];
                std::copy(temp_indices.begin(), temp_indices.begin() + n_this, indices + ix1);
                for (int ix = 0; ix < n_this; ix++)
                    temp_values[ix] = values[argsorted[ix]];
                std::copy(temp_values.begin(), temp_values.begin() + n_this, values + ix1);

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
void sort_sparse_indices_numeric_known_ncol
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    int ncol
)
{
    sort_sparse_indices_known_ncol(
        INTEGER(indptr),
        INTEGER(indices),
        REAL(values),
        indptr.size()-1, ncol
    );
}

// [[Rcpp::export(rng = false)]]
void sort_sparse_indices_logical_known_ncol
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    int ncol
)
{
    sort_sparse_indices_known_ncol(
        INTEGER(indptr),
        INTEGER(indices),
        LOGICAL(values),
        indptr.size()-1, ncol
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

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector deepcopy_num(Rcpp::NumericVector x)
{
    return Rcpp::NumericVector(x.begin(), x.end());
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector deepcopy_int(Rcpp::IntegerVector x)
{
    return Rcpp::IntegerVector(x.begin(), x.end());
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector deepcopy_log(Rcpp::LogicalVector x)
{
    return Rcpp::LogicalVector(x.begin(), x.end());
}

// [[Rcpp::export(rng = false)]]
Rcpp::String deepcopy_str(Rcpp::String x)
{
    return Rcpp::String(x);
}
