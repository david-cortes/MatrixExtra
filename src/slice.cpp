#include "MatrixExtra.h"

size_t get_size_reserve(size_t nnz, size_t take1, size_t take2)
{
    if (sizeof(size_t) < sizeof(uint64_t))
    {
        uint64_t larger = std::min((uint64_t)nnz, (uint64_t)take1 * (uint64_t)take2);
        if (larger >= SIZE_MAX)
            return nnz;
        else
            return larger;
    }

    else
    {
        return std::min(nnz, (size_t)take1 * (size_t)take2);
    }
}

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
bool check_is_rev_seq(Rcpp::IntegerVector indices)
{
    if (indices.size() < 2) return true;
    int n_els = indices.size();
    if ((indices[0] - indices[n_els - 1]) != n_els - 1) return false;
    for (int ix = 1; ix < n_els; ix++) {
        if (indices[ix] != indices[ix - 1] - 1) return false;
    }
    return true;
}

template <class RcppVector, class InputDType>
Rcpp::List reverse_rows_template
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    RcppVector values_
)
{
    Rcpp::IntegerVector indptr_new(indptr.size());
    Rcpp::IntegerVector indices_new_(indices_.size());
    RcppVector values_new_;

    const int *restrict indices = INTEGER(indices_);
    const InputDType *restrict values = nullptr;
    int *restrict indices_new = INTEGER(indices_new_);
    InputDType *restrict values_new = nullptr;
    if (values_.size()) {
        values_new_ = RcppVector(values_.size());
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
            values = (InputDType*)REAL(values_);
            values_new = (InputDType*)REAL(values_new_);
        }
        else if (std::is_same<RcppVector, Rcpp::LogicalVector>::value) {
            values = (InputDType*)INTEGER(values_);
            values_new = (InputDType*)LOGICAL(values_new_);
        }
    }

    int nrows = indptr.size() - 1;
    int rev_row;
    int n_this;
    for (int row = 0; row < nrows; row++)
    {
        rev_row = nrows - row - 1;
        n_this = indptr[rev_row+1] - indptr[rev_row];
        indptr_new[row+1] = indptr_new[row] + n_this;
        std::copy(indices + indptr[rev_row], indices + indptr[rev_row+1], indices_new + indptr_new[row]);
        if (values)
        std::copy(values + indptr[rev_row], values + indptr[rev_row+1], values_new + indptr_new[row]);
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr_new,
        Rcpp::_["indices"] = indices_new_,
        Rcpp::_["values"] = values_new_
    );
}


// [[Rcpp::export(rng = false)]]
Rcpp::List reverse_rows_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values
)
{
    return reverse_rows_template<Rcpp::NumericVector, double>(
        indptr,
        indices,
        values
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List reverse_rows_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values
)
{
    return reverse_rows_template<Rcpp::LogicalVector, int>(
        indptr,
        indices,
        values
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List reverse_rows_binary
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices
)
{
    return reverse_rows_template<Rcpp::NumericVector, double>(
        indptr,
        indices,
        Rcpp::NumericVector()
    );
}

template <class RcppVector>
void reverse_columns_inplace
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    RcppVector values_,
    int ncol
)
{
    int *restrict indices = INTEGER(indices_);
    auto values = values_.begin();
    const bool has_values = values_.size();

    int nrows = indptr.size() - 1;
    for (int row = 0; row < nrows; row++)
    {
        if (indptr[row] < indptr[row+1])
        {
            #pragma omp simd
            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                indices[ix] = ncol - indices[ix] - 1;
            std::reverse(indices + indptr[row], indices + indptr[row+1]);
            if (has_values)
            std::reverse(values + indptr[row], values + indptr[row+1]);
        }
    }
}

// [[Rcpp::export(rng = false)]]
void reverse_columns_inplace_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::NumericVector values_,
    int ncol
)
{
    return reverse_columns_inplace(
        indptr,
        indices_,
        values_,
        ncol
    );
}

// [[Rcpp::export(rng = false)]]
void reverse_columns_inplace_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::LogicalVector values_,
    int ncol
)
{
    return reverse_columns_inplace(
        indptr,
        indices_,
        values_,
        ncol
    );
}

// [[Rcpp::export(rng = false)]]
void reverse_columns_inplace_binary
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::NumericVector values_,
    int ncol
)
{
    return reverse_columns_inplace(
        indptr,
        indices_,
        Rcpp::NumericVector(),
        ncol
    );
}

/* TODO: make these ones work also for logical vectors */

template <class RcppVector>
Rcpp::List copy_csr_rows_template
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppVector values,
    Rcpp::IntegerVector rows_take
)
{
    size_t total_size = 0;
    for (const int row : rows_take) total_size += indptr[row + 1] - indptr[row];
    if (total_size == 0) {
        return Rcpp::List::create(Rcpp::_["indptr"] = Rcpp::IntegerVector(),
                                  Rcpp::_["indices"] = Rcpp::IntegerVector(),
                                  Rcpp::_["values"] = RcppVector());
    }
    
    Rcpp::IntegerVector new_indptr = Rcpp::IntegerVector(rows_take.size() + 1);
    Rcpp::IntegerVector new_indices = Rcpp::IntegerVector(total_size);
    RcppVector new_values = RcppVector(values.size()? total_size : 0);

    size_t n_copy;
    int row;
    int *restrict ptr_indptr = indptr.begin();
    int *restrict ptr_indices = indices.begin();
    auto *restrict prt_values = values.begin();
    int *restrict ptr_new_indptr = new_indptr.begin();
    int *restrict ptr_new_indices = new_indices.begin();
    auto *restrict ptr_new_values = new_values.begin();
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
Rcpp::List copy_csr_rows_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_take
)
{
    return copy_csr_rows_template(
        indptr,
        indices,
        values,
        rows_take
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_rows_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    Rcpp::IntegerVector rows_take
)
{
    return copy_csr_rows_template(
        indptr,
        indices,
        values,
        rows_take
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_rows_binary
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::IntegerVector rows_take
)
{
    return copy_csr_rows_template(
        indptr,
        indices,
        Rcpp::NumericVector(),
        rows_take
    );
}

template <class RcppVector>
Rcpp::List copy_csr_rows_col_seq_template
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppVector values,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take,
    const bool index1
)
{
    const int min_col = *std::min_element(cols_take.begin(), cols_take.end()) - index1;
    const int max_col = *std::max_element(cols_take.begin(), cols_take.end()) - index1;
    Rcpp::IntegerVector new_indptr(rows_take.size() + 1);

    int *restrict ptr_indptr = indptr.begin();
    int *restrict ptr_indices = indices.begin();
    auto *restrict ptr_values = values.begin();
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
    auto *restrict ptr_new_values = new_values.begin();

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
Rcpp::List copy_csr_rows_col_seq_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take,
    const bool index1
)
{
    return copy_csr_rows_col_seq_template(
        indptr,
        indices,
        values,
        rows_take,
        cols_take,
        index1
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_rows_col_seq_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take,
    const bool index1
)
{
    return copy_csr_rows_col_seq_template(
        indptr,
        indices,
        values,
        rows_take,
        cols_take,
        index1
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_rows_col_seq_binary
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take,
    const bool index1
)
{
    return copy_csr_rows_col_seq_template(
        indptr,
        indices,
        Rcpp::NumericVector(),
        rows_take,
        cols_take,
        index1
    );
}

/* TODO: if the indices are sorted, this function could use std::lower_bound instead */

template <class RcppVector, class InputDType, class CompileFlag>
Rcpp::List copy_csr_arbitrary_template
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppVector values,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take
)
{
    Rcpp::IntegerVector new_indptr(rows_take.size()+1);

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

    const int min_j = *std::min_element(cols_take.begin(), cols_take.end());
    const int max_j = *std::max_element(cols_take.begin(), cols_take.end());

    std::vector<int> new_indices;
    std::vector<InputDType> new_values;

    std::vector<int> argsort_cols(cols_take.size());
    std::vector<int> temp_int(cols_take.size());
    std::vector<InputDType> temp_double(values.size()? cols_take.size() : 0);


    size_t size_reserve = get_size_reserve(indices.size(), rows_take.size(), cols_take.size());
    new_indices.reserve(size_reserve);
    if (std::is_same<CompileFlag, bool>::value) new_values.reserve(size_reserve);

    int size_this = 0;
    int row = 0;
    std::unordered_map<int, int>::iterator match;

    for (int row_ix = 0; row_ix < (int)rows_take.size(); row_ix++)
    {
        row = rows_take[row_ix];
        for (int ix = indptr[row]; ix < indptr[row + 1]; ix++)
        {
            if (indices[ix] < min_j || indices[ix] > max_j)
                continue;

            match = new_mapping.find(indices[ix]);
            if (match != new_mapping.end())
            {
                if (has_duplicates && n_repeats[indices[ix]] > 1)
                {
                    for (const auto el : indices_rep[indices[ix]]) {
                        new_indices.push_back(el);
                        if (std::is_same<CompileFlag, bool>::value)
                        new_values.push_back(values[ix]);
                    }
                }
                else
                {
                    new_indices.push_back(match->second);
                    if (std::is_same<CompileFlag, bool>::value)
                    new_values.push_back(values[ix]);
                }
            }
        }
        new_indptr[row_ix + 1] = new_indices.size();
        if (!cols_are_sorted && new_indptr[row_ix + 1] > new_indptr[row_ix])
        {
            size_this = new_indptr[row_ix + 1] - new_indptr[row_ix];
            std::iota(argsort_cols.begin(), argsort_cols.begin() + size_this,
                      new_indptr[row_ix]);
            std::sort(argsort_cols.begin(), argsort_cols.begin() + size_this,
                      [&new_indices](const int a, const int b) {
                        return new_indices[a] < new_indices[b];
                    });
            for (int col = 0; col < size_this; col++) {
                temp_int[col] = new_indices[argsort_cols[col]];
                if (std::is_same<CompileFlag, bool>::value)
                temp_double[col] = new_values[argsort_cols[col]];
            }
            std::copy(temp_int.begin(), temp_int.begin() + size_this,
                      new_indices.begin() + new_indptr[row_ix]);
            if (std::is_same<CompileFlag, bool>::value)
            std::copy(temp_double.begin(), temp_double.begin() + size_this,
                      new_values.begin() + new_indptr[row_ix]);
        }
    }

    Rcpp::List out;
    out["indptr"] = new_indptr;
    VectorConstructorArgs args;
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &new_indices;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    new_indices.clear();
    if (values.size()) {
        if (std::is_same<RcppVector, Rcpp::LogicalVector>::value) {
            args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &new_values;
            args.as_logical = true;
        } else if (std::is_same<RcppVector, Rcpp::IntegerVector>::value) {
            args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &new_values;
        } else {
            args.as_integer = false; args.from_cpp_vec = true; args.num_vec_from = &new_values;
        }
        out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    }
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_arbitrary_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take
)
{
    return copy_csr_arbitrary_template<Rcpp::NumericVector, double, bool>(
        indptr,
        indices,
        values,
        rows_take,
        cols_take
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_arbitrary_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take
)
{
    return copy_csr_arbitrary_template<Rcpp::LogicalVector, int, bool>(
        indptr,
        indices,
        values,
        rows_take,
        cols_take
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List copy_csr_arbitrary_binary
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::IntegerVector rows_take,
    Rcpp::IntegerVector cols_take
)
{
    return copy_csr_arbitrary_template<Rcpp::NumericVector, double, int>(
        indptr,
        indices,
        Rcpp::NumericVector(),
        rows_take,
        cols_take
    );
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

template <class real_t>
real_t extract_single_val_csr
(
    int *restrict indptr,
    int *restrict indices,
    real_t *restrict values,
    const int row, const int col,
    const bool is_sorted
)
{
    if (indptr[row] == indptr[row+1])
        return 0;
    else {
        int *st_this = indices + indptr[row];
        int *end_this = indices + indptr[row+1];
        if (!is_sorted) {
            for (int *ix = st_this; ix < end_this; ix++)
                if (*ix == col)
                    return values? values[ix - indices] : 1;
            return 0;
        }

        if (col < *st_this || col > *(end_this-1)) {
            return 0;
        }

        int *res = std::lower_bound(st_this, end_this, col);
        if (res >= end_this || *res != col) {
            return 0;
        }
        else {
            if (values != nullptr)
                return values[res - indices];
            else
                return 1;
        }
    }

}

double extract_single_val_csr
(
    int *restrict indptr,
    int *restrict indices,
    double *restrict values,
    const int row, const int col,
    const bool is_sorted
)
{
    return extract_single_val_csr<double>(
        indptr,
        indices,
        values,
        row, col,
        is_sorted
    );
}

int extract_single_val_csr
(
    int *restrict indptr,
    int *restrict indices,
    int *restrict values,
    const int row, const int col,
    const bool is_sorted
)
{
    return extract_single_val_csr<int>(
        indptr,
        indices,
        values,
        row, col,
        is_sorted
    );
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
    return extract_single_val_csr<double>(
        INTEGER(indptr),
        INTEGER(indices),
        REAL(values),
        row, col, false
    );
}

// [[Rcpp::export(rng = false)]]
int extract_single_val_csr_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    int row, int col
)
{
    return extract_single_val_csr<int>(
        INTEGER(indptr),
        INTEGER(indices),
        LOGICAL(values),
        row, col, false
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
    return extract_single_val_csr<double>(
        INTEGER(indptr),
        INTEGER(indices),
        (double*)nullptr,
        row, col, false
    );
}
