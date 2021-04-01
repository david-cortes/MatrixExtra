#include "MatrixExtra.h"

// [[Rcpp::export(rng = false)]]
bool check_is_seq(Rcpp::IntegerVector indices)
{
    if (indices.size() < 2) return true;
    int n_els = indices.size();
    if ((indices[n_els - 1] - indices[0]) != n_els - 1) return false;
    for (int ix = 1; ix < n_els; ix++) {
        if (indices[ix] != indices[ix - 1] + 1) return false;
    }
    return true;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_rows(Rcpp::IntegerVector indptr, Rcpp::IntegerVector indices,
                         Rcpp::NumericVector values, Rcpp::IntegerVector rows_take)
{
    size_t total_size = 0;
    for (const int row : rows_take) total_size += indptr[row + 1] - indptr[row];
    if (total_size == 0) {
        return Rcpp::List::create(Rcpp::_["indptr"] = Rcpp::IntegerVector(),
                                  Rcpp::_["indices"] = Rcpp::IntegerVector(),
                                  Rcpp::_["values"] = Rcpp::NumericVector());
    }
    
    Rcpp::IntegerVector new_indptr = Rcpp::IntegerVector(rows_take.size() + 1);
    Rcpp::IntegerVector new_indices = Rcpp::IntegerVector(total_size);
    Rcpp::NumericVector new_values = Rcpp::NumericVector(values.size()? total_size : 0);

    size_t n_copy;
    int row;
    int *restrict ptr_indptr = indptr.begin();
    int *restrict ptr_indices = indices.begin();
    double *restrict prt_values = values.begin();
    int *restrict ptr_new_indptr = new_indptr.begin();
    int *restrict ptr_new_indices = new_indices.begin();
    double *restrict ptr_new_values = new_values.begin();
    const bool has_values = values.size() > 0;

    size_t curr = 0;
    for (int ix = 0; ix < (int)rows_take.size(); ix++)
    {
        row = rows_take[ix];
        n_copy = ptr_indptr[row + 1] - ptr_indptr[row];
        ptr_new_indptr[ix + 1] = ptr_new_indptr[ix] + n_copy;
        if (n_copy) {
            std::copy(ptr_indices + ptr_indptr[row], ptr_indices + ptr_indptr[row + 1],
                      ptr_new_indices + curr);
            if (has_values)
            std::copy(prt_values + ptr_indptr[row], prt_values + ptr_indptr[row + 1],
                      ptr_new_values + curr);
        }
        curr += n_copy;
    }
    return Rcpp::List::create(Rcpp::_["indptr"] = new_indptr,
                              Rcpp::_["indices"] = new_indices,
                              Rcpp::_["values"] = new_values);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_rows_col_seq(Rcpp::IntegerVector indptr, Rcpp::IntegerVector indices,
                                 Rcpp::NumericVector values,
                                 Rcpp::IntegerVector rows_take,
                                 Rcpp::IntegerVector cols_take)
{
    int min_col = *std::min_element(cols_take.begin(), cols_take.end());
    int max_col = *std::max_element(cols_take.begin(), cols_take.end());
    Rcpp::IntegerVector new_indptr(rows_take.size() + 1);

    int *restrict ptr_indptr = indptr.begin();
    int *restrict ptr_indices = indices.begin();
    double *restrict ptr_values = values.begin();
    int *restrict ptr_new_indptr = new_indptr.begin();
    const bool has_values = values.size() > 0;

    size_t total_size = 0;
    for (int row = 0; row < (int)rows_take.size(); row++)
    {
        for (int ix = ptr_indptr[rows_take[row]]; ix < ptr_indptr[rows_take[row] + 1]; ix++)
            total_size += (ptr_indices[ix] >= min_col) && (ptr_indices[ix] <= max_col);
        ptr_new_indptr[row + 1] = total_size;
    }

    if (total_size == 0) {
        return Rcpp::List::create(Rcpp::_["indptr"] = new_indptr,
                                  Rcpp::_["indices"] = Rcpp::IntegerVector(),
                                  Rcpp::_["values"] = Rcpp::NumericVector());
    }

    Rcpp::IntegerVector new_indices = Rcpp::IntegerVector(total_size);
    Rcpp::NumericVector new_values = Rcpp::NumericVector(has_values? total_size : 0);
    int *restrict ptr_new_indices = new_indices.begin();
    double *restrict ptr_new_values = new_values.begin();

    int curr = 0;
    for (int row = 0; row < (int)rows_take.size(); row++)
    {
        for (int ix = ptr_indptr[rows_take[row]]; ix < ptr_indptr[rows_take[row] + 1]; ix++)
        {
            if ((ptr_indices[ix] >= min_col) && (ptr_indices[ix] <= max_col))
            {
                ptr_new_indices[curr] = ptr_indices[ix] - min_col;
                if (has_values)
                ptr_new_values[curr] = ptr_values[ix];
                curr++;
            }
        }
    }
    return Rcpp::List::create(Rcpp::_["indptr"] = new_indptr,
                              Rcpp::_["indices"] = new_indices,
                              Rcpp::_["values"] = new_values);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_arbitrary(Rcpp::IntegerVector indptr, Rcpp::IntegerVector indices,
                              Rcpp::NumericVector values, Rcpp::IntegerVector rows_take,
                              Rcpp::IntegerVector cols_take)
{
    std::unordered_map<int, int> new_mapping;
    for (int col = 0; col < (int)cols_take.size(); col++) new_mapping[cols_take[col]] = col;
    std::unordered_map<int, int> n_repeats;
    for (auto el : cols_take) n_repeats[el]++;
    bool has_duplicates = false;
    for (auto& el : n_repeats) {
        if (el.second > 1) {
            has_duplicates = true;
            break;
        }
    }
    std::unordered_map<int, std::vector<int>> indices_rep;
    if (has_duplicates) {
        for (int col = 0; col < (int)cols_take.size(); col++) {
            if (n_repeats[cols_take[col]] > 1) {
                indices_rep[cols_take[col]].push_back(col);
            }
        }
    }

    bool cols_are_sorted = true;
    for (int ix = 1; ix < (int)cols_take.size(); ix++) {
        if (cols_take[ix] < cols_take[ix - 1]) {
            cols_are_sorted = false;
            break;
        }
    }

    Rcpp::IntegerVector new_indptr;
    VectorConstructorArgs args;
    args.as_integer = true; args.from_cpp_vec = false; args.size = rows_take.size() + 1;
    new_indptr = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    std::vector<int> new_indices;
    std::vector<double> new_values;

    std::vector<int> argsort_cols;
    std::vector<int> temp_int;
    std::vector<double> temp_double;

    const bool has_values = values.size() > 0;

    int size_this = 0;
    int row = 0;
    for (int row_ix = 0; row_ix < (int)rows_take.size(); row_ix++)
    {
        row = rows_take[row_ix];
        for (int ix = indptr[row]; ix < indptr[row + 1]; ix++)
        {
            auto match = new_mapping.find(indices[ix]);
            if (match != new_mapping.end())
            {
                if (has_duplicates && n_repeats[indices[ix]] > 1)
                {
                    for (const auto& el : indices_rep[indices[ix]]) {
                        new_indices.push_back(el);
                        if (has_values)
                        new_values.push_back(values[ix]);
                    }
                }
                else
                {
                    new_indices.push_back(match->second);
                    if (has_values)
                    new_values.push_back(values[ix]);
                }
            }
        }
        new_indptr[row_ix + 1] = new_indices.size();
        if (!cols_are_sorted && new_indptr[row_ix + 1] > new_indptr[row_ix])
        {
            size_this = new_indptr[row_ix + 1] - new_indptr[row_ix];
            if (argsort_cols.size() < (size_t)size_this) {
                argsort_cols.resize(size_this);
                temp_int.resize(size_this);
                temp_double.resize(size_this);
            }
            std::iota(argsort_cols.begin(), argsort_cols.begin() + size_this,
                      new_indptr[row_ix]);
            std::sort(argsort_cols.begin(), argsort_cols.begin() + size_this,
                      [&new_indices](const int a, const int b) {
                        return new_indices[a] < new_indices[b];
                    });
            for (int col = 0; col < size_this; col++) {
                temp_int[col] = new_indices[argsort_cols[col]];
                if (has_values)
                temp_double[col] = new_values[argsort_cols[col]];
            }
            std::copy(temp_int.begin(), temp_int.begin() + size_this,
                      new_indices.begin() + new_indptr[row_ix]);
            if (has_values)
            std::copy(temp_double.begin(), temp_double.begin() + size_this,
                      new_values.begin() + new_indptr[row_ix]);
        }
    }

    Rcpp::List out;
    out["indptr"] = new_indptr;
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &new_indices;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    new_indices.clear();
    if (values.size()) {
        args.as_integer = false; args.from_cpp_vec = true; args.num_vec_from = &new_values;
        out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    }
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector repeat_indices_n_times(Rcpp::IntegerVector indices,
                                           Rcpp::IntegerVector remainder,
                                           int ix_length, int desired_length)
{
    int full_repeats = desired_length / ix_length;
    auto n_indices = indices.size();
    Rcpp::IntegerVector out(n_indices*full_repeats + remainder.size());
    for (int repetition = 0; repetition < full_repeats; repetition++)
    {
        #pragma omp simd
        for (int ix = 0; ix < n_indices; ix++)
            out[ix + n_indices*repetition] = indices[ix] + ix_length*repetition;
    }
    
    #pragma omp simd
    for (int ix = 0; ix < remainder.size(); ix++)
        out[ix + n_indices*full_repeats] = remainder[ix] + ix_length*full_repeats;
    return out;
}

double extract_single_val_csr
(
    int *restrict indptr,
    int *restrict indices,
    double *restrict values,
    int row, int col,
    bool check_sorted
)
{
    if (indptr[row] == indptr[row+1])
        return 0;
    else {
        int st_this = indptr[row];
        size_t n_this = indptr[row+1] - st_this;
        if (check_sorted && !check_is_sorted(indices + st_this, n_this))
        {
            std::vector<size_t> argsorted(n_this);
            std::vector<int> indices_sorted(n_this);
            std::vector<double> values_sorted(n_this);
            int *restrict indices_ = indices + st_this;
            double *restrict values_ = values? (values + st_this) : nullptr;
            std::iota(argsorted.begin(), argsorted.end(), (size_t)0);
            std::sort(argsorted.begin(), argsorted.end(),
                      [&indices_](const size_t a, const size_t b)
                      {return indices_[a] < indices_[b];});
            for (size_t ix = 0; ix < n_this; ix++)
                indices_sorted[ix] = indices_[argsorted[ix]];
            if (values != nullptr)
            for (size_t ix = 0; ix < n_this; ix++)
                values_sorted[ix] = values_[argsorted[ix]];
            std::copy(indices_sorted.begin(), indices_sorted.end(), indices_);
            if (values != nullptr)
            std::copy(values_sorted.begin(), values_sorted.end(), values_);
        }

        int *res = std::lower_bound(indices + st_this, indices + indptr[row+1], col);
        if (res >= indices + indptr[row+1]) {
            return 0;
        }
        else {
            if (values != nullptr)
                return values[res - (indices + st_this)];
            else
                return 1;
        }
    }

}

// [[Rcpp::export(rng = false)]]
double extract_single_val_csr_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    int row, int col
)
{
    return extract_single_val_csr(
        INTEGER(indptr),
        INTEGER(indices),
        REAL(values),
        row, col, true
    );
}

// [[Rcpp::export(rng = false)]]
double extract_single_val_csr_binary
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    int row, int col
)
{
    return extract_single_val_csr(
        INTEGER(indptr),
        INTEGER(indices),
        (double*)nullptr,
        row, col, true
    );
}
