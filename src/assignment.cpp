#include "MatrixExtra.h"

#ifdef __clang__
#   pragma clang diagnostic push
#   pragma clang diagnostic ignored "-Wpass-failed"
#endif

static inline void insert_col_into_row
(
    int *restrict indptr,
    int *restrict indices,
    double *restrict values,
    int *restrict indices_new,
    double *restrict values_new,
    const int row_set,
    const int col_set,
    const double val_set,
    int *restrict argsorted,
    int *restrict buffer,
    int &curr,
    bool &has_col,
    const bool search_for_col
)
{
    has_col = false;
    int pos;

    if (!search_for_col)
        goto insert_col;

    for (pos = indptr[row_set]; pos < indptr[row_set+1]; pos++)
    {
        if (indices[pos] == col_set) {
            has_col = true;
            break;
        }
    }

    if (has_col)
    {
        std::copy(indices + indptr[row_set],
                  indices + indptr[row_set+1],
                  indices_new + curr);
        std::copy(values + indptr[row_set],
                  values + indptr[row_set+1],
                  values_new + curr);
        values_new[curr + (pos - indptr[row_set])] = val_set;
    }

    else
    {
        insert_col:
        if (indptr[row_set+1] == indptr[row_set]) {
            indices_new[curr] = col_set;
            values_new[curr] = val_set;
        }

        else if (col_set < indices[indptr[row_set]]) {
            indices_new[curr] = col_set;
            values_new[curr] = val_set;
            std::copy(indices + indptr[row_set],
                      indices + indptr[row_set+1],
                      indices_new + curr + 1);
            std::copy(values + indptr[row_set],
                      values + indptr[row_set+1],
                      values_new + curr + 1);
            check_and_sort_single_row_inplace(
                indices_new + curr,
                values_new + curr,
                argsorted,
                buffer,
                indptr[row_set+1] - indptr[row_set] + 1,
                true
            );
        }

        else if (col_set > indices[indptr[row_set+1]-1]) {
            std::copy(indices + indptr[row_set],
                      indices + indptr[row_set+1],
                      indices_new + curr);
            std::copy(values + indptr[row_set],
                      values + indptr[row_set+1],
                      values_new + curr);
            indices_new[curr + indptr[row_set+1] - indptr[row_set]] = col_set;
            values_new[curr + indptr[row_set+1] - indptr[row_set]] = val_set;
            check_and_sort_single_row_inplace(
                indices_new + curr,
                values_new + curr,
                argsorted,
                buffer,
                indptr[row_set+1] - indptr[row_set] + 1,
                true
            );
        }

        else {
            indices_new[curr] = col_set;
            values_new[curr] = val_set;
            std::copy(indices + indptr[row_set],
                      indices + indptr[row_set+1],
                      indices_new + curr + 1);
            std::copy(values + indptr[row_set],
                      values + indptr[row_set+1],
                      values_new + curr + 1);
            check_and_sort_single_row_inplace(
                indices_new + curr,
                values_new + curr,
                argsorted,
                buffer,
                indptr[row_set+1] - indptr[row_set] + 1,
                false
            );
        }
    }

    curr += (indptr[row_set+1] - indptr[row_set]) + (int)!has_col;
}

static inline void remove_col_from_row
(
    int *restrict indptr,
    int *restrict indices,
    double *restrict values,
    int *restrict indices_new,
    double *restrict values_new,
    const int row_set,
    const int col_set,
    int &curr,
    bool &has_col
)
{
    has_col = false;
    int pos;
    for (pos = indptr[row_set]; pos < indptr[row_set+1]; pos++)
    {
        if (indices[pos] == col_set)
        {
            has_col = true;
            break;
        }
    }

    if (!has_col)
    {
        std::copy(indices + indptr[row_set], indices + indptr[row_set+1], indices_new + curr);
        std::copy(values + indptr[row_set], values + indptr[row_set+1], values_new + curr);
    }

    else if (pos == indptr[row_set])
    {
        if (indptr[row_set+1] - indptr[row_set] > 1)
        {
            std::copy(indices + indptr[row_set] + 1, indices + indptr[row_set+1], indices_new + curr);
            std::copy(values + indptr[row_set] + 1, values + indptr[row_set+1], values_new + curr);
        }
    }

    else if (pos == indptr[row_set+1])
    {
        if (indptr[row_set+1] - indptr[row_set] > 1)
        {
            std::copy(indices + indptr[row_set], indices + pos, indices_new + curr);
            std::copy(values + indptr[row_set], values + pos, values_new + curr);
        }
    }

    else
    {
        std::copy(indices + indptr[row_set], indices + pos, indices_new + curr);
        std::copy(indices + pos + 1, indices + indptr[row_set+1], indices_new + curr + (pos - indptr[row_set]));

        std::copy(values + indptr[row_set], values + pos, values_new + curr);
        std::copy(values + pos + 1, values + indptr[row_set+1], values_new + curr + (pos - indptr[row_set]));
    }

    curr += (indptr[row_set+1] - indptr[row_set]) - (int)has_col;
}

static inline int sizeof_setdiff
(
    int *restrict s1,
    int *restrict s2,
    int n1, int n2
)
{
    int n_diff = n1;
    int *restrict end1 = s1 + n1;
    int *restrict end2 = s2 + n2;

    while (true)
    {
        if (s1 >= end1 || s2 >= end2) {
            break;
        }

        else if (*s1 == *s2) {
            s1++;
            s2++;
            n_diff--;
        }

        else if (*s1 > *s2) {
            s2 = std::lower_bound(s2, end2, *s1);
        }

        else {
            s1 = std::lower_bound(s1, end1, *s2);
        }
    }

    return n_diff;
}

static inline int sizeof_setintersect
(
    int *restrict s1,
    int *restrict s2,
    int n1, int n2
)
{
    if (n1 == 0 || n2 == 0) return 0;
    int n_intersect = 0;
    int *restrict end1 = s1 + n1;
    int *restrict end2 = s2 + n2;

    while (true)
    {
        if (s1 >= end1 || s2 >= end2) {
            break;
        }

        else if (*s1 == *s2) {
            s1++;
            s2++;
            n_intersect++;
        }

        else if (*s1 > *s2) {
            s2 = std::lower_bound(s2, end2, *s1);
        }

        else {
            s1 = std::lower_bound(s1, end1, *s2);
        }
    }

    return n_intersect;
}

static inline int remove_cols_from_row
(
    int *restrict indices,
    double *restrict values,
    int n,
    int *restrict cols_remove,
    int n_remove,
    int *restrict indices_new,
    double *restrict values_new
)
{
    if (n == 0) return 0;
    int *restrict end_indices = indices + n;
    int *restrict end_cols = cols_remove + n_remove;
    int *temp, diff;
    int n_removed = 0;

    while (true)
    {
        if (indices >= end_indices || cols_remove >= end_cols) {
            std::copy(indices, end_indices, indices_new);
            std::copy(values, values + (end_indices - indices), values_new);
            break;
        }

        else if (*indices == *cols_remove) {
            indices++;
            values++;
            cols_remove++;
            n_removed++;
        }

        else if (*indices > *cols_remove) {
            cols_remove = std::lower_bound(cols_remove, end_cols, *indices);
        }

        else {
            temp = std::lower_bound(indices, end_indices, *cols_remove);
            diff = temp - indices;
            std::copy(indices, temp, indices_new);
            std::copy(values, values + diff, values_new);
            indices = temp;
            values += diff;
            indices_new += diff;
            values_new += diff;
        }
    }

    return n_removed;
}

static inline int set_cols_in_row_to_const_value
(
    int *restrict indices,
    double *restrict values,
    int n,
    int *restrict cols_add,
    int n_add,
    const double val_set,
    int *restrict indices_new,
    double *restrict values_new
)
{
    int *restrict end_indices = indices + n;
    int *restrict end_cols = cols_add + n_add;
    int *temp, diff;
    int n_tot = 0;

    if (n == 0)
    {
        std::copy(cols_add, cols_add + n_add, indices_new);
        std::fill_n(values_new, n_add, val_set);
        return n_add;
    }

    while (true)
    {
        if (indices >= end_indices || cols_add >= end_cols) {
            if (indices < end_indices) {
                std::copy(indices, end_indices, indices_new);
                std::copy(values, values + (end_indices - indices), values_new);
                n_tot += end_indices - indices;
            }
            for (; cols_add < end_cols; cols_add++) {
                *indices_new = *cols_add;
                *values_new = val_set;
                indices_new++;
                values_new++;
                n_tot++;
            }
            break;
        }

        else if (*indices == *cols_add) {
            *indices_new = *cols_add;
            *values_new = val_set;
            indices++;
            values++;
            cols_add++;
            indices_new++;
            values_new++;
            n_tot++;
        }

        else if (*indices > *cols_add) {
            for (; cols_add < end_cols && *cols_add < *indices; cols_add++) {
                *indices_new = *cols_add;
                *values_new = val_set;
                indices_new++;
                values_new++;
                n_tot++;
            }
        }

        else {
            temp = std::lower_bound(indices, end_indices, *cols_add);
            diff = temp - indices;
            std::copy(indices, temp, indices_new);
            std::copy(values, values + diff, values_new);
            indices = temp;
            values += diff;
            indices_new += diff;
            values_new += diff;
            n_tot += diff;
        }
    }

    return n_tot;
}

#define check_max_size(diff, indices) \
    if ((diff) >= INT_MAX - (indices).size()) \
        Rcpp::stop("Error: resulting matrix would be larger than INT_MAX limit.\n")

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_row_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int row_set
)
{
    int n_this = indptr[row_set+1] - indptr[row_set];
    if (n_this == 0) {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size() - n_this);
    Rcpp::NumericVector values_new(indices.size() - n_this);

    const int nrows = indptr.size()-1;

    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (int row = row_set; row < nrows; row++)
        indptr_new[row+1] -= n_this;

    std::copy(indices.begin(), indices.begin() + indptr[row_set], indices_new.begin());
    std::copy(indices.begin() + indptr[row_set+1], indices.end(), indices_new.begin() + indptr[row_set]);
    std::copy(values.begin(), values.begin() + indptr[row_set], values_new.begin());
    std::copy(values.begin() + indptr[row_set+1], values.end(), values_new.begin() + indptr[row_set]);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_col_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int col_set
)
{
    int n_affected = 0;
    for (auto el : indices)
        n_affected += el == col_set;

    if (n_affected == 0) {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size() - n_affected);
    Rcpp::NumericVector values_new(indices.size() - n_affected);

    int n_this;
    int cum_affected = 0;
    int curr = 0;

    const int nrows = indptr.size() - 1;

    for (int row = 0; row < nrows; row++)
    {
        n_this = 0;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            n_this += indices[ix] == col_set;

        if (!n_this)
        {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[row+1],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[row+1],
                      values_new.begin() + indptr_new[row]);
            curr += indptr[row+1] - indptr[row];
        }

        else
        {
            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            {
                if (indices[ix] != col_set)
                {
                    indices_new[curr] = indices[ix];
                    values_new[curr] = values[ix];
                    curr++;
                }
            }
        }

        cum_affected += n_this;
        indptr_new[row+1] -= cum_affected;
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_row_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int ncols,
    const int row_set,
    const double val_set
)
{
    const int n_before = indptr[row_set+1] - indptr[row_set];
    const int diff = ncols - n_before;

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());
        std::fill_n(values_new.begin() + indptr[row_set], ncols, val_set);

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    const int nrows = indptr.size()-1;

    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (int row = row_set; row < nrows; row++)
        indptr_new[row+1] += diff;

    std::copy(indices.begin(),
              indices.begin() + indptr[row_set],
              indices_new.begin());
    std::iota(indices_new.begin() + indptr[row_set],
              indices_new.begin() + indptr[row_set] + ncols,
              0);
    std::copy(indices.begin() + indptr[row_set+1],
              indices.end(),
              indices_new.begin() + indptr[row_set] + ncols);

    std::copy(values.begin(),
              values.begin() + indptr[row_set],
              values_new.begin());
    std::fill_n(values_new.begin() + indptr[row_set], ncols, val_set);
    std::copy(values.begin() + indptr[row_set+1],
              values.end(),
              values_new.begin() + indptr[row_set] + ncols);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_col_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int ncols,
    const int col_set,
    const double val_set
)
{
    const int nrows = indptr.size()-1;
    int n_prev = 0;
    for (auto el : indices)
        n_prev += el == col_set;
    const int diff = nrows - n_prev;

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());
        const int nnz = indices.size();
        for (int ix = 0; ix < nnz; ix++)
            values_new[ix] = (indices[ix] == col_set)? val_set : values[ix];

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);

    int curr = 0;
    int cum_affected = 0;
    bool has_col;

    for (int row = 0; row < nrows; row++)
    {
        insert_col_into_row(
            indptr.begin(),
            indices.begin(),
            values.begin(),
            indices_new.begin(),
            values_new.begin(),
            row,
            col_set,
            val_set,
            argsorted.get(),
            buffer.get(),
            curr,
            has_col,
            true
        );
        cum_affected += (int)!has_col;

        indptr_new[row+1] += cum_affected;
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_val_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int row_set,
    const int col_set
)
{
    bool has_val = false;
    int ix;
    for (ix = indptr[row_set]; ix < indptr[row_set+1]; ix++)
    {
        if (indices[ix] == col_set)
        {
            has_val = true;
            break;
        }
    }

    if (!has_val)
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size()-1);
    Rcpp::NumericVector values_new(indices.size()-1);

    const int nrows = indptr.size()-1;

    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (int row = row_set; row < nrows; row++)
        indptr_new[row+1] -= 1;

    std::copy(indices.begin(), indices.begin() + ix, indices_new.begin());
    std::copy(indices.begin() + ix + 1, indices.end(), indices_new.begin() + ix);
    std::copy(values.begin(), values.begin() + ix, values_new.begin());
    std::copy(values.begin() + ix + 1, values.end(), values_new.begin() + ix);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_val_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int ncols,
    const int row_set,
    const int col_set,
    const double val_set
)
{
    bool has_val = false;
    int ix;
    for (ix = indptr[row_set]; ix < indptr[row_set+1]; ix++)
    {
        if (indices[ix] == col_set)
        {
            has_val = true;
            break;
        }
    }

    if (has_val)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());
        values_new[ix] = val_set;
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size()+1);
    Rcpp::NumericVector values_new(indices.size()+1);

    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);

    const int nrows = indptr.size()-1;

    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (int row = row_set; row < nrows; row++)
        indptr_new[row+1] += 1;

    std::copy(indices.begin(), indices.begin() + indptr[row_set], indices_new.begin());
    std::copy(values.begin(), values.begin() + indptr[row_set], values_new.begin());

    int curr = indptr[row_set];
    bool unused;

    insert_col_into_row(
        indptr.begin(),
        indices.begin(),
        values.begin(),
        indices_new.begin(),
        values_new.begin(),
        row_set,
        col_set,
        val_set,
        argsorted.get(),
        buffer.get(),
        curr,
        unused,
        false
    );

    std::copy(indices.begin() + indptr[row_set+1],
              indices.end(),
              indices_new.begin() + indptr[row_set+1] + 1);
    std::copy(values.begin() + indptr[row_set+1],
              values.end(),
              values_new.begin() + indptr[row_set+1] + 1);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_row_to_rowvec
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int ncols,
    const int row_set,
    Rcpp::NumericVector vec_set
)
{
    const int diff = ncols - (indptr[row_set+1] - indptr[row_set]);
    const int n_repeats = ncols / vec_set.size();

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());
        for (int repetition = 0; repetition < n_repeats; repetition++)
            std::copy(vec_set.begin(), vec_set.end(),
                      values_new.begin() + indptr[row_set] + repetition*vec_set.size());

        if (check_is_sorted(indices.begin() + indptr[row_set], ncols)) {
            return Rcpp::List::create(
                Rcpp::_["indptr"] = indptr,
                Rcpp::_["indices"] = indices,
                Rcpp::_["values"] = values_new
            );
        }

        else {
            Rcpp::IntegerVector indices_new(indices.begin(), indices.end());
            std::iota(indices_new.begin() + indptr[row_set], indices_new.begin() + indptr[row_set] + ncols, 0);
            return Rcpp::List::create(
                Rcpp::_["indptr"] = indptr,
                Rcpp::_["indices"] = indices_new,
                Rcpp::_["values"] = values_new
            );
        }
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    const int nrows = indptr.size() - 1;

    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (int row = row_set; row < nrows; row++)
        indptr_new[row+1] += diff;

    std::copy(indices.begin(), indices.begin() + indptr[row_set], indices_new.begin());
    std::iota(indices_new.begin() + indptr[row_set], indices_new.begin() + indptr[row_set] + ncols, 0);
    std::copy(indices.begin() + indptr[row_set+1], indices.end(), indices_new.begin() + indptr_new[row_set+1]);

    std::copy(values.begin(), values.begin() + indptr[row_set], values_new.begin());
    for (int repetition = 0; repetition < n_repeats; repetition++)
        std::copy(vec_set.begin(), vec_set.end(), values_new.begin() + indptr[row_set] + repetition*vec_set.size());
    std::copy(values.begin() + indptr[row_set+1], values.end(), values_new.begin() + indptr_new[row_set+1]);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_col_to_colvec
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int ncols,
    const int col_set,
    Rcpp::NumericVector vec_set
)
{
    const int nrows = indptr.size() - 1;
    int n_prev = 0;
    for (auto el : indices)
        n_prev += el == col_set;
    const int diff = nrows - n_prev;

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());
        for (int row = 0; row < nrows; row++)
        {
            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            {
                if (indices[ix] == col_set)
                    values_new[ix] = vec_set[row % vec_set.size()];
            }
        }

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);

    bool has_col;
    int curr = 0;
    int cum_affected = 0;

    for (int row = 0; row < nrows; row++)
    {
        insert_col_into_row(
            indptr.begin(),
            indices.begin(),
            values.begin(),
            indices_new.begin(),
            values_new.begin(),
            row,
            col_set,
            vec_set[row % vec_set.size()],
            argsorted.get(),
            buffer.get(),
            curr,
            has_col,
            true
        );
        cum_affected += (int)!has_col;

        indptr_new[row+1] += cum_affected;
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_row_to_svec
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int ncols,
    const int row_set,
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const int length
)
{
    if (indices.size() == 0 && ii.size() == 0)
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    const int n_repeats = ncols / length;
    const int nnz = ii.size();
    if (((indptr[row_set+1] - indptr[row_set]) % nnz) == 0 && indptr[row_set] < indptr[row_set+1])
    {

        if (n_repeats <= 1)
        {
            for (int ix = indptr[row_set]; ix < indptr[row_set+1]; ix++)
                if (indices[ix] != ii[ix - indptr[row_set]])
                    goto normal_route;
        }

        else
        {
            for (int ix = indptr[row_set]; ix < indptr[row_set+1]; ix++)
                if (indices[ix] != (ii[(ix - indptr[row_set]) % nnz] + length * (indices[ix] / length)))
                    goto normal_route;
        }

        /* If reaching here, it's the same indices */
        Rcpp::NumericVector values_new(values.begin(), values.end());
        for (int repetition = 0; repetition < n_repeats; repetition++)
            std::copy(xx.begin(), xx.end(), values_new.begin() + indptr[row_set] + repetition*nnz);

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    normal_route:
    const int diff = nnz*n_repeats - (indptr[row_set+1] - indptr[row_set]);
    Rcpp::IntegerVector indptr_new(indptr.begin(), indptr.end());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(values.size() + diff);

    const int nrows = indptr.size() - 1;

    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (int row = row_set; row < nrows; row++)
        indptr_new[row+1] += diff;

    std::copy(indices.begin(), indices.begin() + indptr[row_set], indices_new.begin());
    for (int repetition = 0; repetition < n_repeats; repetition++) {
        std::copy(ii.begin(), ii.end(), indices_new.begin() + indptr[row_set] + nnz*repetition);
        if (repetition > 0) {
            #ifdef _OPENMP
            #pragma omp simd
            #endif
            for (int ix = 0; ix < nnz; ix++)
                indices_new[indptr[row_set] + nnz*repetition + ix] += length*repetition;
        }
    }
    if (row_set < nrows-1)
        std::copy(indices.begin() + indptr[row_set+1], indices.end(), indices_new.begin() + indptr[row_set+1] + diff);

    std::copy(values.begin(), values.begin() + indptr[row_set], values_new.begin());
    for (int repetition = 0; repetition < n_repeats; repetition++)
        std::copy(xx.begin(), xx.end(), values_new.begin() + indptr[row_set] + nnz*repetition);
    if (row_set < nrows-1)
        std::copy(values.begin() + indptr[row_set+1], values.end(), values_new.begin() + indptr[row_set+1] + diff);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_col_to_svec
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int ncols,
    const int col_set,
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const int length
)
{
    const int nrows = indptr.size() - 1;
    const int n_repeats = nrows / length;
    const int nnz = ii.size();

    Rcpp::IntegerVector indptr_new(nrows+1);
    std::unique_ptr<int[]> indices_new(new int[indices.size() + nnz*n_repeats]);
    std::unique_ptr<double[]> values_new(new double[indices.size() + nnz*n_repeats]);

    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);
    

    int curr = 0;
    int row = 0;
    const auto begin_i = ii.begin();
    const auto end_i = ii.end();
    auto curr_i = begin_i;
    int row_end;
    int offset;
    bool has_col;

    for (int repetition = 0; repetition < n_repeats; repetition++)
    {
        offset = repetition * length;
        row = offset;
        row_end = offset + length;
        curr_i = begin_i;

        while (true)
        {
            if (curr_i >= end_i || row >= row_end) {

                for (; row < row_end; row++)
                {
                    remove_col_from_row(
                        indptr.begin(),
                        indices.begin(),
                        values.begin(),
                        indices_new.get(),
                        values_new.get(),
                        row,
                        col_set,
                        curr,
                        has_col
                    );
                    indptr_new[row+1] = indptr[row+1] - indptr[row] - (int)has_col;
                }
                break;

            }

            else if (row == *curr_i + offset) {

                insert_col_into_row(
                    indptr.begin(),
                    indices.begin(),
                    values.begin(),
                    indices_new.get(),
                    values_new.get(),
                    row,
                    col_set,
                    xx[curr_i - begin_i],
                    argsorted.get(),
                    buffer.get(),
                    curr,
                    has_col,
                    true
                );
                indptr_new[row+1] = indptr[row+1] - indptr[row] + (int)!has_col;

                row++;
                curr_i++;

            }

            else {

                for (; row < *curr_i + offset; row++)
                {
                    remove_col_from_row(
                        indptr.begin(),
                        indices.begin(),
                        values.begin(),
                        indices_new.get(),
                        values_new.get(),
                        row,
                        col_set,
                        curr,
                        has_col
                    );
                    indptr_new[row+1] = indptr[row+1] - indptr[row] - (int)has_col;
                }

            }
        }
    }

    for (int row = 0; row < nrows; row++)
        indptr_new[row+1] += indptr_new[row];


    Rcpp::List out;
    out["indptr"] = indptr_new;
    VectorConstructorArgs args;
    args.from_pointer = true; args.as_integer = true; args.size = curr;
    args.int_pointer_from = indices_new.get();
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_new.reset();
    args.as_integer = false; args.num_pointer_from = values_new.get();
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_rowseq_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int row_set_st,
    const int row_set_end
)
{
    const int nrows = indptr.size() - 1;
    const int diff = indptr[row_set_end+1] - indptr[row_set_st];
    const int n_set = row_set_end - row_set_st + 1;
    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() - diff);
    Rcpp::NumericVector values_new(indices.size() - diff);

    std::copy(indices.begin(), indices.begin() + indptr[row_set_st], indices_new.begin());
    std::copy(indices.begin() + indptr[row_set_end+1], indices.end(), indices_new.begin() + indptr[row_set_st]);

    std::copy(values.begin(), values.begin() + indptr[row_set_st], values_new.begin());
    std::copy(values.begin() + indptr[row_set_end+1], values.end(), values_new.begin() + indptr[row_set_st]);

    std::copy(indptr.begin(), indptr.begin() + row_set_st + 1, indptr_new.begin());
    std::fill_n(indptr_new.begin() + row_set_st + 1, n_set, indptr[row_set_st]);
    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (int row = row_set_end; row < nrows; row++)
        indptr_new[row+1] = indptr[row+1] - diff;

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_rowseq_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int row_set_st,
    const int row_set_end,
    const int ncols,
    const double val_set
)
{
    const int nrows = indptr.size() - 1;
    const int n_fill = ncols * (row_set_end - row_set_st + 1);
    const int diff = n_fill - (indptr[row_set_end+1] - indptr[row_set_st]);
    check_max_size(diff, indices);

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());
        std::fill(values_new.begin() + indptr[row_set_st], values_new.begin() + indptr[row_set_end+1], val_set);
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    Rcpp::IntegerVector indptr_new(nrows+1);
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    std::copy(indices.begin(), indices.begin() + indptr[row_set_st], indices_new.begin());
    int curr = indptr[row_set_st];
    for (int row = row_set_st; row <= row_set_end; row++) {
        std::iota(indices_new.begin() + curr, indices_new.begin() + curr + ncols, 0);
        curr += ncols;
    }
    std::copy(indices.begin() + indptr[row_set_end+1], indices.end(), indices_new.begin() + curr);

    std::copy(values.begin(), values.begin() + indptr[row_set_st], values_new.begin());
    std::fill(values_new.begin() + indptr[row_set_st],
              values_new.begin() + indptr[row_set_st] + n_fill,
              val_set);
    std::copy(values.begin() + indptr[row_set_end+1],
              values.end(),
              values_new.begin() + indptr[row_set_st] + n_fill);

    std::copy(indptr.begin(), indptr.begin() + row_set_st + 1, indptr_new.begin());
    for (int row = row_set_st; row <= row_set_end; row++) {
        indptr_new[row+1] = indptr_new[row] + ncols;
    }
    for (int row = row_set_end+1; row < nrows; row++)
        indptr_new[row+1] = indptr_new[row] + indptr[row+1] - indptr[row];

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_colseq_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int col_set_st,
    const int col_set_end,
    const int ncols
)
{
    const int nrows = indptr.size() - 1;
    const int ncurr = std::accumulate(indices.begin(), indices.end(), 0,
                                      [&col_set_st, &col_set_end](const int tot, const int el)
                                      {return tot + (el >= col_set_st && el <= col_set_end);});
    if (ncurr == 0) {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }


    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() - ncurr);
    Rcpp::NumericVector values_new(indices.size() - ncurr);


    std::unique_ptr<int[]> ix_take(new int[ncols]);
    int curr = 0;
    int n_this;
    for (int row = 0; row < nrows; row++)
    {
        n_this = 0;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
        {
            if (indices[ix] < col_set_st || indices[ix] > col_set_end)
                ix_take[n_this++] = ix;
        }
        if (n_this)
        {
            for (int ix = 0; ix < n_this; ix++)
                indices_new[curr + ix] = indices[ix_take[ix]] ;
            for (int ix = 0; ix < n_this; ix++)
                values_new[curr + ix] = values[ix_take[ix]] ;
            curr += n_this;
        }
        indptr_new[row+1] = curr;
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_colseq_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int col_set_st,
    const int col_set_end,
    const int ncols,
    const double val_set
)
{
    const int ncurr = std::accumulate(indices.begin(), indices.end(), 0,
                                      [&col_set_st, &col_set_end](const int tot, const int el)
                                      {return tot + (el >= col_set_st && el <= col_set_end);});
    const int nrows = indptr.size() - 1;
    const int n_fill = col_set_end - col_set_st + 1;
    const int diff = nrows * n_fill - ncurr;
    check_max_size(diff, indices);

    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    std::unique_ptr<int[]> ix_take(new int[ncols]);
    int curr = 0;
    int n_left, n_right;
    for (int row = 0; row < nrows; row++)
    {
        n_left = 0;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
        {
            if (indices[ix] < col_set_st)
                ix_take[n_left++] = ix;
        }
        if (n_left)
        {
            for (int ix = 0; ix < n_left; ix++)
                indices_new[curr + ix] = indices[ix_take[ix]] ;
            for (int ix = 0; ix < n_left; ix++)
                values_new[curr + ix] = values[ix_take[ix]] ;
            curr += n_left;
        }

        std::iota(indices_new.begin() + curr, indices_new.begin() + curr + n_fill, col_set_st);
        std::fill(values_new.begin() + curr, values_new.begin() + curr + n_fill, val_set);
        curr += n_fill;

        n_right = 0;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
        {
            if (indices[ix] > col_set_end)
                ix_take[n_right++] = ix;
        }
        if (n_right)
        {
            for (int ix = 0; ix < n_right; ix++)
                indices_new[curr + ix] = indices[ix_take[ix]] ;
            for (int ix = 0; ix < n_right; ix++)
                values_new[curr + ix] = values[ix_take[ix]] ;
            curr += n_right;
        }

        indptr_new[row+1] = curr;
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_rows_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_set
)
{
    const int nrows = indptr.size() - 1;
    std::sort(rows_set.begin(), rows_set.end());
    const int diff = std::accumulate(rows_set.begin(), rows_set.end(), 0,
                                     [&indptr](const int tot, const int row)
                                     {return tot + (indptr[row+1] - indptr[row]);});
    if (diff == 0) {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() - diff);
    Rcpp::NumericVector values_new(indices.size() - diff);

    int row = 0;
    auto curr_ix = rows_set.begin();
    const auto end_ix = rows_set.end();
    int curr = 0;
    while (true)
    {
        if (curr_ix >= end_ix || row >= nrows) {
            if (row < nrows)
            {
                std::copy(indices.begin() + indptr[row],
                          indices.begin() + indptr[nrows],
                          indices_new.begin() + curr);
                std::copy(values.begin() + indptr[row],
                          values.begin() + indptr[nrows],
                          values_new.begin() + curr);
                for (; row < nrows; row++)
                    indptr_new[row+1] = indptr[row+1] - indptr[row];
            }
            break;
        }

        else if (row == *curr_ix) {
            row++;
            curr_ix++;
        }

        else if (row > *curr_ix) {
            curr_ix = std::lower_bound(curr_ix, end_ix, row);
        }

        else {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[*curr_ix],
                      indices_new.begin() + curr);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[*curr_ix],
                      values_new.begin() + curr);
            curr += indptr[*curr_ix] - indptr[row];
            for (; row < *curr_ix; row++)
                indptr_new[row+1] = indptr[row+1] - indptr[row];
        }
    }

    for (int row = 0; row < nrows; row++)
        indptr_new[row+1] += indptr_new[row];

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_rows_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_set,
    const int ncols,
    const double val_set
)
{
    const int nrows = indptr.size() - 1;
    std::sort(rows_set.begin(), rows_set.end());
    const int ncurr = std::accumulate(rows_set.begin(), rows_set.end(), 0,
                                      [&indptr](const int tot, const int row)
                                      {return tot + (indptr[row+1] - indptr[row]);});
    const int diff = rows_set.size()*ncols - ncurr;
    

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());
        for (int row : rows_set)
            std::fill(values_new.begin() + indptr[row], values_new.begin() + indptr[row+1], val_set);

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    
    else if (ncurr == 0)
    {
        Rcpp::IntegerVector indptr_new(indptr.size());
        Rcpp::IntegerVector indices_new(indices.size() + diff);
        Rcpp::NumericVector values_new(indices.size() + diff);

        auto curr_ix = rows_set.begin();
        int last_ix = rows_set[rows_set.size()-1];
        for (int row = 0; row <= last_ix; row++)
        {
            if (row == *curr_ix)
            {
                std::iota(indices_new.begin() + indptr_new[row],
                          indices_new.begin() + indptr_new[row] + ncols,
                          0);
                std::fill_n(values_new.begin() + indptr_new[row], ncols, val_set);
                indptr_new[row+1] = indptr_new[row] + ncols;
                curr_ix++;
            }

            else
            {
                std::copy(indices.begin() + indptr[row],
                          indices.begin() + indptr[row+1],
                          indices_new.begin() + indptr_new[row]);
                std::copy(values.begin() + indptr[row],
                          values.begin() + indptr[row+1],
                          values_new.begin() + indptr_new[row]);
                indptr_new[row+1] = indptr_new[row] + indptr[row+1] - indptr[row];

            }
        }

        if (last_ix < nrows - 1)
        {
            std::copy(indices.begin() + indptr[last_ix+1],
                      indices.end(),
                      indices_new.begin() + indptr_new[last_ix+1]);
            std::copy(values.begin() + indptr[last_ix+1],
                      values.end(),
                      values_new.begin() + indptr_new[last_ix+1]);
            for (int row = last_ix+1; row < nrows; row++)
                indptr_new[row+1] = indptr_new[row] + indptr[row+1] - indptr[row];
        }

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr_new,
            Rcpp::_["indices"] = indices_new,
            Rcpp::_["values"] = values_new
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);


    int row = 0;
    auto curr_ix = rows_set.begin();
    const auto end_ix = rows_set.end();
    int curr = 0;
    while (true)
    {
        if (curr_ix >= end_ix || row >= nrows) {
            if (row < nrows)
            {
                std::copy(indices.begin() + indptr[row],
                          indices.begin() + indptr[nrows],
                          indices_new.begin() + curr);
                std::copy(values.begin() + indptr[row],
                          values.begin() + indptr[nrows],
                          values_new.begin() + curr);
                for (; row < nrows; row++)
                    indptr_new[row+1] = indptr[row+1] - indptr[row];
            }
            break;
        }

        else if (row == *curr_ix) {
            std::iota(indices_new.begin() + curr, indices_new.begin() + curr + ncols, 0);
            std::fill_n(values_new.begin() + curr, ncols, val_set);
            indptr_new[row+1] = ncols;

            curr += ncols;
            row++;
            curr_ix++;
        }

        else if (row > *curr_ix) {
            for (; *curr_ix < row; curr_ix++) {
                std::iota(indices_new.begin() + curr, indices_new.begin() + curr + ncols, 0);
                std::fill_n(values_new.begin() + curr, ncols, val_set);
                indptr_new[*curr_ix+1] = ncols;
                curr += ncols;
            }
        }

        else {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[*curr_ix],
                      indices_new.begin() + curr);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[*curr_ix],
                      values_new.begin() + curr);
            curr += indptr[*curr_ix] - indptr[row];
            for (; row < *curr_ix; row++)
                indptr_new[row+1] = indptr[row+1] - indptr[row];
        }
    }


    for (int row = 0; row < nrows; row++)
        indptr_new[row+1] += indptr_new[row];

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_cols_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector cols_set,
    const int ncols
)
{
    const int nrows = indptr.size() - 1;
    hashed_set<int> cols_set_(cols_set.begin(), cols_set.end());
    const int col_min = *std::min_element(cols_set.begin(), cols_set.end());
    const int col_max = *std::max_element(cols_set.begin(), cols_set.end());

    const int ncurr = std::accumulate(indices.begin(), indices.end(), 0,
                                      [&cols_set_, &col_min, &col_max](const int tot, const int el)
                                      {return tot + (el >= col_min &&
                                                     el <= col_max &&
                                                     cols_set_.find(el) != cols_set_.end());});
    if (ncurr == 0)
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    VectorConstructorArgs args;
    args.as_integer = true; args.size = indptr.size();
    Rcpp::IntegerVector indptr_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    args.size = indices.size() - ncurr;
    Rcpp::IntegerVector indices_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    args.as_integer = false;
    Rcpp::NumericVector values_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);

    std::unique_ptr<int[]> ix_take(new int[ncols]);
    int n_this;
    for (int row = 0; row < nrows; row++)
    {
        n_this = 0;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
        {
            if (indices[ix] < col_min || indices[ix] > col_max || cols_set_.find(indices[ix]) == cols_set_.end())
                ix_take[n_this++] = ix;
        }
        if (n_this)
        {
            for (int ix = 0; ix < n_this; ix++)
                indices_new[indptr_new[row] + ix] = indices[ix_take[ix]];
            for (int ix = 0; ix < n_this; ix++)
                values_new[indptr_new[row] + ix] = values[ix_take[ix]];
        }
        indptr_new[row+1] = indptr_new[row] + n_this;
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_cols_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector cols_set,
    const int ncols,
    const double val_set
)
{
    const int nrows = indptr.size() - 1;
    std::sort(cols_set.begin(), cols_set.end());
    hashed_set<int> cols_set_(cols_set.begin(), cols_set.end());
    const int col_min = *std::min_element(cols_set.begin(), cols_set.end());
    const int col_max = *std::max_element(cols_set.begin(), cols_set.end());

    const int ncurr = std::accumulate(indices.begin(), indices.end(), 0,
                                      [&cols_set_, &col_min, &col_max](const int tot, const int el)
                                      {return tot + (el >= col_min &&
                                                     el <= col_max &&
                                                     cols_set_.find(el) != cols_set_.end());});
    const int diff = nrows*cols_set.size() - ncurr;
    const int n_fill = cols_set.size();

    VectorConstructorArgs args;

    if (diff == 0)
    {
        args.as_integer = false; args.size = values.size(); args.from_pointer = true;
        args.num_pointer_from = (void*)((double*)(values.begin()));
        Rcpp::NumericVector values_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);

        for (int row = 0; row < nrows; row++)
        {
            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
            {
                if (indices[ix] >= col_min && indices[ix] <= col_max &&
                    cols_set_.find(indices[ix]) != cols_set_.end())
                {
                    values_new[ix] = val_set;
                }
            }
        }

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    check_max_size(diff, indices);

    args.as_integer = true; args.size = indptr.size();
    Rcpp::IntegerVector indptr_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    args.size = indices.size() + diff;
    Rcpp::IntegerVector indices_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    args.as_integer = false;
    Rcpp::NumericVector values_new = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);

    std::unique_ptr<int[]> ix_take(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);

    if (ncurr == 0)
    {
        for (int row = 0; row < nrows; row++)
        {
            std::copy(cols_set.begin(),
                      cols_set.end(),
                      indices_new.begin() + indptr_new[row]);
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[row+1],
                      indices_new.begin() + indptr_new[row] + n_fill);
            std::fill_n(values_new.begin() + indptr_new[row], n_fill, val_set);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[row+1],
                      values_new.begin() + indptr_new[row] + n_fill);
            check_and_sort_single_row_inplace(
                indices_new.begin() + indptr_new[row],
                values_new.begin() + indptr_new[row],
                ix_take.get(),
                buffer.get(),
                n_fill + indptr[row+1] - indptr[row],
                indptr[row] == indptr[row+1]
            );
            indptr_new[row+1] = indptr_new[row] + n_fill + indptr[row+1] - indptr[row];
        }

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr_new,
            Rcpp::_["indices"] = indices_new,
            Rcpp::_["values"] = values_new
        );
    }

    int n_this;
    
    for (int row = 0; row < nrows; row++)
    {
        n_this = 0;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
        {
            if (indices[ix] < col_min || indices[ix] > col_max || cols_set_.find(indices[ix]) == cols_set_.end())
                ix_take[n_this++] = ix;
        }
        for (int ix = 0; ix < n_this; ix++)
            indices_new[indptr_new[row] + ix] = indices[ix_take[ix]];
        std::copy(cols_set.begin(), cols_set.end(), indices_new.begin() + indptr_new[row] + n_this);
        for (int ix = 0; ix < n_this; ix++)
            values_new[indptr_new[row] + ix] = values[ix_take[ix]];
        std::fill_n(values_new.begin() + indptr_new[row] + n_this, n_fill, val_set);

        check_and_sort_single_row_inplace(
            indices_new.begin() + indptr_new[row],
            values_new.begin() + indptr_new[row],
            ix_take.get(),
            buffer.get(),
            n_fill + n_this,
            indptr[row] == indptr[row+1]
        );

        indptr_new[row+1] = indptr_new[row] + n_fill + n_this;
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_rows_single_col_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_set,
    const int col_set,
    const int ncols
)
{
    Rcpp::IntegerVector indptr_new(indptr.size());
    int ncurr = 0;
    for (int row : rows_set)
    {
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
        {
            ncurr += indices[ix] == col_set;
            indptr_new[row+1] -= indices[ix] == col_set;
        }
    }

    if (ncurr == 0)
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    const int nrows = indptr.size() - 1;
    for (int row = 0; row < nrows; row++)
        indptr_new[row+1] += indptr_new[row] + indptr[row+1] - indptr[row];

    Rcpp::IntegerVector indices_new(indices.size() - ncurr);
    Rcpp::NumericVector values_new(indices.size() - ncurr);
    std::sort(rows_set.begin(), rows_set.end());

    int row = 0;
    auto curr_ix = rows_set.begin();
    const auto end_ix = rows_set.end();

    std::unique_ptr<int[]> ix_take(new int[ncols]);
    int n_this;

    while (true)
    {
        if (curr_ix >= end_ix || row >= nrows) {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[nrows],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[nrows],
                      values_new.begin() + indptr_new[row]);
            break;
        }

        else if (row == *curr_ix) {
            if (indptr[row] < indptr[row+1])
            {
                if ((indptr_new[row+1] - indptr_new[row]) == (indptr[row+1] - indptr[row])) {
                    std::copy(indices.begin() + indptr[row],
                              indices.begin() + indptr[row+1],
                              indices_new.begin() + indptr_new[row]);
                    std::copy(values.begin() + indptr[row],
                              values.begin() + indptr[row+1],
                              values_new.begin() + indptr_new[row]);
                }

                else {
                    n_this = 0;
                    for (int ix = indptr[row]; ix < indptr[row+1]; ix++) {
                        if (indices[ix] != col_set)
                            ix_take[n_this++] = ix;
                    }
                    for (int ix = 0; ix < n_this; ix++) {
                        indices_new[ix + indptr_new[row]] = indices[ix_take[ix]];
                        values_new[ix + indptr_new[row]] = values[ix_take[ix]];
                    }
                }
            }

            row++;
            curr_ix++;
        }

        else {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[*curr_ix],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[*curr_ix],
                      values_new.begin() + indptr_new[row]);
            row = *curr_ix;
        }
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_rows_single_col_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_set,
    const int col_set,
    const double val_set,
    const int ncols
)
{
    const int nrows = indptr.size() - 1;
    Rcpp::IntegerVector indptr_new(indptr.size());
    int diff = 0;
    bool has_col;
    for (int row : rows_set)
    {
        has_col = false;
        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
        {
            if (indices[ix] == col_set)
            {
                has_col = true;
                break;
            }
        }
        indptr_new[row+1] = !has_col;
        diff += !has_col;
    }

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());

        for (int row : rows_set)
            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                values_new[ix] = (indices[ix] == col_set)? val_set : values_new[ix];

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);

    std::sort(rows_set.begin(), rows_set.end());

    for (int row = 0; row < nrows; row++)
        indptr_new[row+1] += indptr_new[row] + indptr[row+1] - indptr[row];

    int row = 0;
    auto curr_ix = rows_set.begin();
    const auto end_ix = rows_set.end();

    while (true)
    {
        if (curr_ix >= end_ix || row >= nrows) {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[nrows],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[nrows],
                      values_new.begin() + indptr_new[row]);
            break;
        }

        else if (row == *curr_ix) {
            if (indptr[row] < indptr[row+1])
            {
                if ((indptr_new[row+1] - indptr_new[row]) == (indptr[row+1] - indptr[row])) {
                    std::copy(indices.begin() + indptr[row],
                              indices.begin() + indptr[row+1],
                              indices_new.begin() + indptr_new[row]);
                    for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                        values_new[indptr_new[row] - indptr[row] + ix]
                            =
                        (indices[ix] == col_set)? val_set : values[ix];
                }

                else {

                    if (col_set < indices[indptr[row]])
                    {
                        indices_new[indptr_new[row]] = col_set;
                        std::copy(indices.begin() + indptr[row],
                                  indices.begin() + indptr[row+1],
                                  indices_new.begin() + indptr_new[row] + 1);
                        values_new[indptr_new[row]] = val_set;
                        std::copy(values.begin() + indptr[row],
                                  values.begin() + indptr[row+1],
                                  values_new.begin() + indptr_new[row] + 1);
                    }

                    else if (col_set > indices[indptr[row+1]-1])
                    {
                        std::copy(indices.begin() + indptr[row],
                                  indices.begin() + indptr[row+1],
                                  indices_new.begin() + indptr_new[row]);
                        indices_new[indptr_new[row+1]-1] = col_set;
                        std::copy(values.begin() + indptr[row],
                                  values.begin() + indptr[row+1],
                                  values_new.begin() + indptr_new[row]);
                        values_new[indptr_new[row+1]-1] = val_set;
                    }

                    else
                    {
                        indices_new[indptr_new[row]] = col_set;
                        std::copy(indices.begin() + indptr[row],
                                  indices.begin() + indptr[row+1],
                                  indices_new.begin() + indptr_new[row] + 1);
                        values_new[indptr_new[row]] = val_set;
                        std::copy(values.begin() + indptr[row],
                                  values.begin() + indptr[row+1],
                                  values_new.begin() + indptr_new[row] + 1);
                        check_and_sort_single_row_inplace(
                            indices_new.begin() + indptr_new[row],
                            values_new.begin() + indptr_new[row],
                            argsorted.get(),
                            buffer.get(),
                            indptr_new[row+1] - indptr_new[row],
                            false
                        );
                    }
                }
            }

            else
            {
                indices_new[indptr_new[row]] = col_set;
                values_new[indptr_new[row]] = val_set;
            }

            row++;
            curr_ix++;
        }

        else {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[*curr_ix],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[*curr_ix],
                      values_new.begin() + indptr_new[row]);
            row = *curr_ix;
        }
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_row_arbitrary_cols_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int row_set,
    Rcpp::IntegerVector cols_set,
    const int ncols
)
{
    const int n_this = indptr[row_set+1] - indptr[row_set];
    if (n_this == 0)
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);
    std::sort(cols_set.begin(), cols_set.end());
    check_and_sort_single_row_inplace(
        indices.begin() + indptr[row_set],
        values.begin() + indptr[row_set],
        argsorted.get(),
        buffer.get(),
        n_this,
        true
    );
    argsorted.reset();
    buffer.reset();

    const int diff = sizeof_setdiff(
        indices.begin() + indptr[row_set],
        cols_set.begin(),
        n_this,
        cols_set.size()
    );

    if (diff == n_this)
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    else if (diff == 0)
    {
        return set_single_row_to_zero(
            indptr,
            indices,
            values,
            row_set
        );
    }

    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() - n_this + diff);
    Rcpp::NumericVector values_new(indices.size() - n_this + diff);
    const int nrows = indptr.size()-1;

    std::copy(indptr.begin(), indptr.begin() + row_set + 1, indptr_new.begin());
    indptr_new[row_set+1] = indptr_new[row_set] + diff;
    for (int row = row_set+1; row < nrows; row++)
        indptr_new[row+1] = indptr_new[row] + indptr[row+1] - indptr[row];

    std::copy(indices.begin(), indices.begin() + indptr[row_set], indices_new.begin());
    remove_cols_from_row(
        indices.begin() + indptr[row_set],
        values.begin() + indptr[row_set],
        n_this,
        cols_set.begin(),
        cols_set.size(),
        indices_new.begin() + indptr_new[row_set],
        values_new.begin() + indptr_new[row_set]
    );
    if (row_set < nrows-1)
        std::copy(indices.begin() + indptr[row_set+1],
                  indices.end(),
                  indices_new.begin() + indptr_new[row_set+1]);

    std::copy(values.begin(), values.begin() + indptr[row_set], values_new.begin());
    if (row_set < nrows-1)
        std::copy(values.begin() + indptr[row_set+1],
                  values.end(),
                  values_new.begin() + indptr_new[row_set+1]);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_single_row_arbitrary_cols_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int row_set,
    Rcpp::IntegerVector cols_set,
    const int ncols,
    const double val_set
)
{
    std::sort(cols_set.begin(), cols_set.end());
    const int size_before = indptr[row_set+1] - indptr[row_set];
    const int n_common = sizeof_setintersect(
        indices.begin() + indptr[row_set],
        cols_set.begin(),
        size_before, cols_set.size()
    );
    const int new_size = ((int)cols_set.size() - n_common) + size_before;

    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);
    check_and_sort_single_row_inplace(
        indices.begin() + indptr[row_set],
        values.begin() + indptr[row_set],
        argsorted.get(),
        buffer.get(),
        size_before,
        true
    );
    argsorted.reset();
    buffer.reset();

    if (new_size == size_before)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());

        set_cols_in_row_to_const_value(
            indices.begin() + indptr[row_set],
            values.begin() + indptr[row_set],
            size_before,
            cols_set.begin(),
            cols_set.size(),
            val_set,
            indices.begin() + indptr[row_set],
            values_new.begin() + indptr[row_set]
        );

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    const int nrows = indptr.size()-1;
    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() + (new_size - size_before));
    Rcpp::NumericVector values_new(indices.size() + (new_size - size_before));

    std::copy(indptr.begin(), indptr.begin() + row_set + 1, indptr_new.begin());
    indptr_new[row_set+1] = indptr_new[row_set] + new_size;
    for (int row = row_set+1; row < nrows; row++)
        indptr_new[row+1] = indptr_new[row] + indptr[row+1] - indptr[row];

    std::copy(indices.begin(), indices.begin() + indptr[row_set], indices_new.begin());
    set_cols_in_row_to_const_value(
        indices.begin() + indptr[row_set],
        values.begin() + indptr[row_set],
        size_before,
        cols_set.begin(),
        cols_set.size(),
        val_set,
        indices_new.begin() + indptr[row_set],
        values_new.begin() + indptr[row_set]
    );
    if (row_set < nrows-1)
        std::copy(indices.begin() + indptr[row_set+1], indices.end(), indices_new.begin() + indptr_new[row_set+1]);
    
    std::copy(values.begin(), values.begin() + indptr[row_set], values_new.begin());
    if (row_set < nrows-1)
        std::copy(values.begin() + indptr[row_set+1], values.end(), values_new.begin() + indptr_new[row_set+1]);


    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_rows_arbitrary_cols_to_zero
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_set,
    Rcpp::IntegerVector cols_set,
    const int ncols
)
{
    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);
    Rcpp::IntegerVector indptr_new(indptr.size());
    std::sort(cols_set.begin(), cols_set.end());

    int ncurr = 0;
    int n_this;
    for (int row : rows_set)
    {
        n_this = indptr[row+1] - indptr[row];
        if (n_this)
        {
            check_and_sort_single_row_inplace(
                indices.begin() + indptr[row],
                values.begin() + indptr[row],
                argsorted.get(),
                buffer.get(),
                n_this,
                true
            );
            indptr_new[row+1] = sizeof_setdiff(
                indices.begin() + indptr[row],
                cols_set.begin(),
                n_this,
                cols_set.size()
            ) - n_this;
            ncurr += -indptr_new[row+1];
        }
    }

    if (ncurr == 0)
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    argsorted.reset();
    buffer.reset();

    const int nrows = indptr.size()-1;
    Rcpp::IntegerVector indices_new(indices.size() - ncurr);
    Rcpp::NumericVector values_new(indices.size() - ncurr);
    std::sort(rows_set.begin(), rows_set.end());

    for (int row = 0; row < nrows; row++)
        indptr_new[row+1] += indptr_new[row] + indptr[row+1] - indptr[row];

    int row = 0;
    auto curr_ix = rows_set.begin();
    const auto end_ix = rows_set.end();
    while (true)
    {
        if (curr_ix >= end_ix || row >= nrows) {
            if (row < nrows) {
                std::copy(indices.begin() + indptr[row], indices.end(), indices_new.begin() + indptr_new[row]);
                std::copy(values.begin() + indptr[row], values.end(), values_new.begin() + indptr_new[row]);
            }
            break;
        }

        else if (*curr_ix == row) {
            remove_cols_from_row(
                indices.begin() + indptr[row],
                values.begin() + indptr[row],
                indptr[row+1] - indptr[row],
                cols_set.begin(),
                cols_set.size(),
                indices_new.begin() + indptr_new[row],
                values_new.begin() + indptr_new[row]
            );
            curr_ix++;
            row++;
        }

        else {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[*curr_ix],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[*curr_ix],
                      values_new.begin() + indptr_new[row]);
            row = *curr_ix;
        }
    }


    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_rows_arbitrary_cols_to_const
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_set,
    Rcpp::IntegerVector cols_set,
    const int ncols,
    const double val_set
)
{
    std::unique_ptr<int[]> argsorted(new int[ncols]);
    std::unique_ptr<int[]> buffer(new int[size_times_ratio_dbl(ncols)]);
    Rcpp::IntegerVector indptr_new(indptr.size());
    std::sort(cols_set.begin(), cols_set.end());
    const int n_fill = cols_set.size();

    int diff = 0;
    int n_this;
    for (int row : rows_set)
    {
        n_this = indptr[row+1] - indptr[row];
        check_and_sort_single_row_inplace(
            indices.begin() + indptr[row],
            values.begin() + indptr[row],
            argsorted.get(),
            buffer.get(),
            n_this,
            true
        );
        indptr_new[row+1] = -sizeof_setintersect(
            indices.begin() + indptr[row],
            cols_set.begin(),
            n_this, n_fill
        ) + n_fill;
        diff += indptr_new[row+1];
    }

    if (diff == 0)
    {
        Rcpp::NumericVector values_new(values.begin(), values.end());

        for (int row : rows_set)
        {
            set_cols_in_row_to_const_value(
                indices.begin() + indptr[row],
                values.begin() + indptr[row],
                indptr[row+1] - indptr[row],
                cols_set.begin(),
                n_fill,
                val_set,
                indices.begin() + indptr[row],
                values_new.begin() + indptr[row]
            );
        }

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values_new
        );
    }

    const int nrows = indptr.size()-1;
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);
    std::sort(rows_set.begin(), rows_set.end());

    for (int row = 0; row < nrows; row++)
        indptr_new[row+1] += indptr_new[row] + indptr[row+1] - indptr[row];

    int row = 0;
    auto curr_ix = rows_set.begin();
    const auto end_ix = rows_set.end();
    while (true)
    {
        if (curr_ix >= end_ix || row >= nrows) {
            if (row < nrows) {
                std::copy(indices.begin() + indptr[row], indices.end(), indices_new.begin() + indptr_new[row]);
                std::copy(values.begin() + indptr[row], values.end(), values_new.begin() + indptr_new[row]);
            }
            break;
        }

        else if (*curr_ix == row) {
            set_cols_in_row_to_const_value(
                indices.begin() + indptr[row],
                values.begin() + indptr[row],
                indptr[row+1] - indptr[row],
                cols_set.begin(),
                n_fill,
                val_set,
                indices_new.begin() + indptr_new[row],
                values_new.begin() + indptr_new[row]
            );
            curr_ix++;
            row++;
        }

        else {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[*curr_ix],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[*curr_ix],
                      values_new.begin() + indptr_new[row]);
            row = *curr_ix;
        }
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_rowseq_to_smat
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    const int row_set_st,
    const int row_set_end,
    Rcpp::IntegerVector indptr_other,
    Rcpp::IntegerVector indices_other,
    Rcpp::NumericVector values_other
)
{
    const int n_prev = indptr[row_set_end+1] - indptr[row_set_st];
    const int n_new = indptr_other[indptr_other.size()-1];
    const int diff = n_new - n_prev;
    const int nrows = indptr.size()-1;

    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    std::copy(indptr.begin(), indptr.begin() + row_set_st + 1, indptr_new.begin());
    for (int row = row_set_st; row <= row_set_end; row++)
        indptr_new[row+1] = indptr_new[row] + (indptr_other[row-row_set_st+1] - indptr_other[row-row_set_st]);
    for (int row = row_set_end+1; row < nrows; row++)
        indptr_new[row+1] = indptr_new[row] + (indptr[row+1] - indptr[row]);

    std::copy(indices.begin(), indices.begin() + indptr[row_set_st], indices_new.begin());
    std::copy(indices_other.begin(), indices_other.end(), indices_new.begin() + indptr[row_set_st]);
    if (row_set_end < nrows-1)
        std::copy(indices.begin() + indptr[row_set_end+1], indices.end(), indices_new.begin() + indptr_new[row_set_end+1]);

    std::copy(values.begin(), values.begin() + indptr[row_set_st], values_new.begin());
    std::copy(values_other.begin(), values_other.end(), values_new.begin() + indptr[row_set_st]);
    if (row_set_end < nrows-1)
        std::copy(values.begin() + indptr[row_set_end+1], values.end(), values_new.begin() + indptr_new[row_set_end+1]);

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List set_arbitrary_rows_to_smat
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_set,
    Rcpp::IntegerVector indptr_other,
    Rcpp::IntegerVector indices_other,
    Rcpp::NumericVector values_other
)
{
    int ncurr = 0;
    for (auto row : rows_set)
        ncurr += indptr[row+1] - indptr[row];
    const int diff = indptr_other[indptr_other.size()-1] - ncurr;
    const int nrows = indptr.size()-1;

    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new(indices.size() + diff);
    Rcpp::NumericVector values_new(indices.size() + diff);

    const auto begin_ix = rows_set.begin();
    auto curr_ix = begin_ix;
    const auto end_ix = rows_set.end();
    int row = 0;
    int ix_this;

    while (true)
    {
        if (curr_ix >= end_ix || row >= nrows) {
            if (row < nrows-1)
            {
                std::copy(indices.begin() + indptr[row],
                          indices.end(),
                          indices_new.begin() + indptr_new[row]);
                std::copy(values.begin() + indptr[row],
                          values.end(),
                          values_new.begin() + indptr_new[row]);
                for (; row < nrows; row++)
                    indptr_new[row+1] = indptr_new[row] + (indptr[row+1] - indptr[row]);
            }
            break;
        }

        else if (row == *curr_ix) {
            ix_this = curr_ix - begin_ix;
            indptr_new[row+1] = indptr_new[row] + (indptr_other[ix_this + 1] - indptr_other[ix_this]);
            std::copy(indices_other.begin() + indptr_other[ix_this],
                      indices_other.begin() + indptr_other[ix_this+1],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values_other.begin() + indptr_other[ix_this],
                      values_other.begin() + indptr_other[ix_this+1],
                      values_new.begin() + indptr_new[row]);
            row++;
            curr_ix++;
        }

        else {
            std::copy(indices.begin() + indptr[row],
                      indices.begin() + indptr[*curr_ix],
                      indices_new.begin() + indptr_new[row]);
            std::copy(values.begin() + indptr[row],
                      values.begin() + indptr[*curr_ix],
                      values_new.begin() + indptr_new[row]);
            for (; row < *curr_ix; row++)
                indptr_new[row+1] = indptr_new[row] + (indptr[row+1] - indptr[row]);
        }
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new,
        Rcpp::_["values"] = values_new
    );
}

#ifdef __clang__
#   pragma clang diagnostic pop
#endif
