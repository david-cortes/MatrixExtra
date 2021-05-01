#include "MatrixExtra.h"

SEXP SafeRcppVector(void *args_)
{
    VectorConstructorArgs *args = (VectorConstructorArgs*)args_;
    std::vector<int> *int_vec_from = (std::vector<int>*)args->int_vec_from;
    std::vector<double> *num_vec_from = (std::vector<double>*)args->num_vec_from;
    
    if (args->as_integer)
    {
        if (args->from_cpp_vec) {
            if (!args->as_logical) {
                if (!args->cpp_lim_size)
                    return Rcpp::IntegerVector(int_vec_from->begin(), int_vec_from->end());
                else
                    return Rcpp::IntegerVector(int_vec_from->begin(), int_vec_from->begin() + args->size);
            }
            else {
                if (!args->cpp_lim_size)
                    return Rcpp::LogicalVector(int_vec_from->begin(), int_vec_from->end());
                else
                    return Rcpp::LogicalVector(int_vec_from->begin(), int_vec_from->begin() + args->size);
            }
        }

        else if (args->from_pointer) {
            int *int_pointer_from = (int*)args->int_pointer_from;
            if (!args->as_logical) {
                return Rcpp::IntegerVector(int_pointer_from, int_pointer_from + args->size);
            }
            else {
                return Rcpp::LogicalVector(int_pointer_from, int_pointer_from + args->size);
            }
        }

        else {
            if (!args->as_logical)
                return Rcpp::IntegerVector(args->size);
            else
                return Rcpp::LogicalVector(args->size);
        }
    }

    else
    {
        if (args->from_cpp_vec) {
            if (!args->cpp_lim_size)
                return Rcpp::NumericVector(num_vec_from->begin(), num_vec_from->end());
            else
                return Rcpp::NumericVector(num_vec_from->begin(), num_vec_from->begin() + args->size);
        }

        else if (args->from_pointer) {
            double *num_pointer_from = (double*)args->num_pointer_from;
            return Rcpp::NumericVector(num_pointer_from, num_pointer_from + args->size);
        }

        else {
            return Rcpp::NumericVector(args->size);
        }
    }
}

// [[Rcpp::export(rng = false)]]
bool contains_any_zero(Rcpp::NumericVector x)
{
    for (auto el : x)
        if (el == 0)
            return true;
    return false;
}

// [[Rcpp::export(rng = false)]]
bool contains_any_inf(Rcpp::NumericVector x)
{
    for (auto el : x)
        if (isinf(el))
            return true;
    return false;
}

// [[Rcpp::export(rng = false)]]
bool contains_any_neg(Rcpp::NumericVector x)
{
    for (auto el : x)
        if (el < 0)
            return true;
    return false;
}

bool contains_any_nas_or_inf(Rcpp::NumericVector x)
{
    for (auto el : x)
        if (ISNAN(el) || isinf(el))
            return true;
    return false;
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

// [[Rcpp::export(rng = false)]]
bool check_is_sorted(Rcpp::IntegerVector x)
{
    return check_is_sorted(x.begin(), x.size());
}

void check_and_sort_single_row_inplace
(
    int *restrict indices,
    double *restrict values,
    int *restrict argsorted,
    int *restrict buffer,
    const int n,
    const bool pre_check
)
{
    if (!pre_check || !check_is_sorted(indices, n))
    {
        std::iota(argsorted, argsorted + n, 0);
        std::sort(argsorted, argsorted + n,
                  [&indices](const int a, const int b)
                  {return indices[a] < indices[b];});
        for (int ix = 0; ix < n; ix++)
            buffer[ix] = indices[argsorted[ix]];
        memcpy(indices, buffer, (size_t)n*sizeof(int));
        double *buffer_ = (double*)buffer;
        for (int ix = 0; ix < n; ix++)
            buffer_[ix] = values[argsorted[ix]];
        memcpy(values, buffer_, (size_t)n*sizeof(double));
    }
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
                if ((int)argsorted.size() < n_this) {
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

template <class RcppVector, class InputDType>
Rcpp::List remove_zero_valued_csr
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppVector values,
    const bool remove_NAs
)
{
    if (!remove_NAs)
    {
        for (auto el : values)
            if (!el)
                goto remove_zeros;
    }

    else
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            for (auto el : values)
                if (!el || ISNAN(el))
                    goto remove_zeros;
        }

        else
        {
            for (auto el : values)
                if (!el || el == NA_LOGICAL)
                    goto remove_zeros;
        }
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr,
        Rcpp::_["indices"] = indices,
        Rcpp::_["values"] = values
    );

    remove_zeros:
    Rcpp::IntegerVector indptr_new(indptr.size());
    std::unique_ptr<int[]> indices_new(new int[indices.size()]);
    std::unique_ptr<InputDType[]> values_new(new InputDType[values.size()]);
    int nrows = indptr.size() - 1;

    int curr = 0;
    if (!remove_NAs)
    {
        for (int row = 0; row < nrows; row++)
        {
            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            {
                if (values[ix])
                {
                    indices_new[curr] = indices[ix];
                    values_new[curr] = values[ix];
                    curr++;
                }
            }
            indptr_new[row+1] = curr;
        }
    }

    else
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            for (int row = 0; row < nrows; row++)
            {
                for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                {
                    if (values[ix] && !ISNAN(values[ix]))
                    {
                        indices_new[curr] = indices[ix];
                        values_new[curr] = values[ix];
                        curr++;
                    }
                }
                indptr_new[row+1] = curr;
            }
        }

        else
        {
            for (int row = 0; row < nrows; row++)
            {
                for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                {
                    if (values[ix] != NA_LOGICAL)
                    {
                        indices_new[curr] = indices[ix];
                        values_new[curr] = values[ix];
                        curr++;
                    }
                }
                indptr_new[row+1] = curr;
            }
        }
    }

    Rcpp::List out;
    out["indptr"] = indptr_new;
    VectorConstructorArgs args;
    args.as_integer = true; args.from_pointer = true;
    args.size = curr; args.int_pointer_from = indices_new.get();
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_new.reset();
    args.as_integer = false; args.from_pointer = true;
    args.num_pointer_from = values_new.get();
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}


// [[Rcpp::export(rng = false)]]
Rcpp::List remove_zero_valued_csr_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const bool remove_NAs
)
{
    return remove_zero_valued_csr<Rcpp::NumericVector, double>(
        indptr,
        indices,
        values,
        remove_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List remove_zero_valued_csr_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    const bool remove_NAs
)
{
    return remove_zero_valued_csr<Rcpp::LogicalVector, int>(
        indptr,
        indices,
        values,
        remove_NAs
    );
}

template <class RcppVector, class InputDType>
Rcpp::List remove_zero_valued_coo
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    RcppVector xx,
    const bool remove_NAs
)
{
    size_t nnz = ii.size();
    if (!remove_NAs)
    {
        for (auto el : xx)
            if (!el)
                goto remove_zeros;
    }

    else
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            for (auto el : xx)
                if (!el || ISNAN(el))
                    goto remove_zeros;
        }

        else
        {
            for (auto el : xx)
                if (!el || el == NA_LOGICAL)
                    goto remove_zeros;
        }
    }

    return Rcpp::List::create(
        Rcpp::_["ii"] = ii,
        Rcpp::_["jj"] = jj,
        Rcpp::_["xx"] = xx
    );

    remove_zeros:
    std::unique_ptr<size_t[]> take(new size_t[nnz]);
    size_t curr = 0;
    if (!remove_NAs)
    {
        for (size_t ix = 0; ix < nnz; ix++)
            if (xx[ix])
                take[curr++] = ix;
    }

    else
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            for (size_t ix = 0; ix < nnz; ix++)
                if (xx[ix] && !ISNAN(xx[ix]))
                    take[curr++] = ix;
        }

        else
        {
            for (size_t ix = 0; ix < nnz; ix++)
                if (xx[ix] && xx[ix] != NA_LOGICAL)
                    take[curr++] = ix;
        }
    }

    VectorConstructorArgs args;
    args.as_integer = true; args.size = curr;
    Rcpp::IntegerVector ii_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    Rcpp::IntegerVector jj_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
        args.as_integer = false;
    } else {
        args.as_integer = true; args.as_logical = true;
    }
    RcppVector xx_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);

    for (size_t ix = 0; ix < curr; ix++) ii_new[ix] = ii[take[ix]];
    for (size_t ix = 0; ix < curr; ix++) jj_new[ix] = jj[take[ix]];
    for (size_t ix = 0; ix < curr; ix++) xx_new[ix] = xx[take[ix]];

    return Rcpp::List::create(
        Rcpp::_["ii"] = ii_new,
        Rcpp::_["jj"] = jj_new,
        Rcpp::_["xx"] = xx_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List remove_zero_valued_coo_numeric
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::NumericVector xx,
    const bool remove_NAs
)
{
    return remove_zero_valued_coo<Rcpp::NumericVector, double>(
        ii,
        jj,
        xx,
        remove_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List remove_zero_valued_coo_logical
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::LogicalVector xx,
    const bool remove_NAs
)
{
    return remove_zero_valued_coo<Rcpp::LogicalVector, int>(
        ii,
        jj,
        xx,
        remove_NAs
    );
}

template <class RcppVector, class InputDType>
Rcpp::List remove_zero_valued_svec
(
    Rcpp::IntegerVector ii,
    RcppVector xx,
    const bool remove_NAs
)
{
    size_t nnz = ii.size();
    if (!remove_NAs)
    {
        for (auto el : xx)
            if (!el)
                goto remove_zeros;
    }

    else
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            for (auto el : xx)
                if (!el || ISNAN(el))
                    goto remove_zeros;
        }

        else if (std::is_same<RcppVector, Rcpp::LogicalVector>::value)
        {
            for (auto el : xx)
                if (!el || el == NA_LOGICAL)
                    goto remove_zeros;
        }

        else
        {
            for (auto el : xx)
                if (!el || el == NA_INTEGER)
                    goto remove_zeros;
        }
    }

    return Rcpp::List::create(
        Rcpp::_["ii"] = ii,
        Rcpp::_["xx"] = xx
    );

    remove_zeros:
    std::unique_ptr<size_t[]> take(new size_t[nnz]);
    size_t curr = 0;

    if (!remove_NAs)
    {
        for (size_t ix = 0; ix < nnz; ix++)
            if (xx[ix])
                take[curr++] = ix;
    }

    else
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            for (size_t ix = 0; ix < nnz; ix++)
                if (xx[ix])
                    take[curr++] = ix;
        }

        else if (std::is_same<RcppVector, Rcpp::LogicalVector>::value)
        {
            for (size_t ix = 0; ix < nnz; ix++)
                if (xx[ix] && xx[ix] != NA_LOGICAL)
                    take[curr++] = ix;
        }

        else
        {
            for (size_t ix = 0; ix < nnz; ix++)
                if (xx[ix] && xx[ix] != NA_INTEGER)
                    take[curr++] = ix;
        }
    }

    VectorConstructorArgs args;
    args.as_integer = true; args.size = curr;
    Rcpp::IntegerVector ii_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
        args.as_integer = false;
    } else if (std::is_same<RcppVector, Rcpp::LogicalVector>::value) {
        args.as_integer = true; args.as_logical = true;
    } else {
        args.as_integer = true;
    }
    Rcpp::IntegerVector xx_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);

    for (size_t ix = 0; ix < curr; ix++) ii_new[ix] = ii[take[ix]];
    for (size_t ix = 0; ix < curr; ix++) xx_new[ix] = xx[take[ix]];

    return Rcpp::List::create(
        Rcpp::_["ii"] = ii_new,
        Rcpp::_["xx"] = xx_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List remove_zero_valued_svec_numeric
(
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const bool remove_NAs
)
{
    return remove_zero_valued_svec<Rcpp::NumericVector, double>(
        ii,
        xx,
        remove_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List remove_zero_valued_svec_integer
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector xx,
    const bool remove_NAs
)
{
    return remove_zero_valued_svec<Rcpp::IntegerVector, int>(
        ii,
        xx,
        remove_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List remove_zero_valued_svec_logical
(
    Rcpp::IntegerVector ii,
    Rcpp::LogicalVector xx,
    const bool remove_NAs
)
{
    return remove_zero_valued_svec<Rcpp::LogicalVector, int>(
        ii,
        xx,
        remove_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List check_valid_csr_matrix
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    int nrows, int ncols
)
{
    int imin = *std::min_element(indices.begin(), indices.end());
    if (imin < 0) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has negative indices.")
        );
    }
    int imax = *std::max_element(indices.begin(), indices.end());
    if (imax >= ncols) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has invalid column indices.")
        );
    }

    for (auto el : indices) {
        if (el == NA_INTEGER) {
            return Rcpp::List::create(
                Rcpp::_["err"] = Rcpp::String("Matrix has indices with missing values.")
            );
        }
    }

    for (auto el : indptr) {
        if (el == NA_INTEGER) {
            return Rcpp::List::create(
                Rcpp::_["err"] = Rcpp::String("Matrix has missing values in the index pointer.")
            );
        }
    }

    for (int ix = 0; ix < nrows; ix++) {
        if (indptr[ix] > indptr[ix+1]) {
            return Rcpp::List::create(
                Rcpp::_["err"] = Rcpp::String("Matrix index pointer is not monotonicaly increasing.")
            );
        }
    }

    return Rcpp::List();
}

// [[Rcpp::export(rng = false)]]
Rcpp::List check_valid_coo_matrix
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    int nrows, int ncols
)
{
    int imin = *std::min_element(ii.begin(), ii.end());
    if (imin < 0) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has negative indices.")
        );
    }
    int imax = *std::max_element(ii.begin(), ii.end());
    if (imax >= nrows) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has invalid column indices.")
        );
    }
    for (auto el : ii) {
        if (el == NA_INTEGER) {
            return Rcpp::List::create(
                Rcpp::_["err"] = Rcpp::String("Matrix has indices with missing values.")
            );
        }
    }

    int jmin = *std::min_element(jj.begin(), jj.end());
    if (jmin < 0) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has negative indices.")
        );
    }
    int jmax = *std::max_element(jj.begin(), jj.end());
    if (jmax >= ncols) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has invalid column indices.")
        );
    }
    for (auto el : jj) {
        if (el == NA_INTEGER) {
            return Rcpp::List::create(
                Rcpp::_["err"] = Rcpp::String("Matrix has indices with missing values.")
            );
        }
    }

    return Rcpp::List();
}

// [[Rcpp::export(rng = false)]]
Rcpp::List check_valid_svec
(
    Rcpp::IntegerVector ii,
    int nrows
)
{
    int imin = *std::min_element(ii.begin(), ii.end());
    if (imin < 0) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has negative indices.")
        );
    }
    int imax = *std::max_element(ii.begin(), ii.end());
    if (imax >= nrows) {
        return Rcpp::List::create(
            Rcpp::_["err"] = Rcpp::String("Matrix has invalid column indices.")
        );
    }
    for (auto el : ii) {
        if (el == NA_INTEGER) {
            return Rcpp::List::create(
                Rcpp::_["err"] = Rcpp::String("Matrix has indices with missing values.")
            );
        }
    }

    return Rcpp::List(); 
}
