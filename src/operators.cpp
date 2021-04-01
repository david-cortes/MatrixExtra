#include "MatrixExtra.h"

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csr_elemwise(Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
                                 Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
                                 Rcpp::NumericVector values1, Rcpp::NumericVector values2)
{
    if (indptr1.size() == indptr2.size() &&
        indices1.size() == indices2.size() &&
        INTEGER(indptr1) == INTEGER(indptr2) &&
        INTEGER(indices1) == INTEGER(indices2)
    ) {
        Rcpp::NumericVector values_out(values1.size());
        #pragma omp simd
        for (int el = 0; el < (int)values1.size(); el++)
            values_out[el] = values1[el] * values2[el];
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr1,
            Rcpp::_["indices"] = indices1,
            Rcpp::_["values"] = values_out
        );
    }

    size_t nrows = indptr1.size() - 1;
    Rcpp::List out = Rcpp::List::create(
        Rcpp::_["indptr"] = Rcpp::IntegerVector(nrows+1)
    );

    int *indptr_out = INTEGER(out["indptr"]);
    std::vector<int> indices_out;
    std::vector<double> values_out;
    indices_out.reserve(std::min(indices1.size(), indices2.size()));
    values_out.reserve(std::min(indices1.size(), indices2.size()));

    int *ptr_indices1 = INTEGER(indices1);
    int *ptr_indices2 = INTEGER(indices2);

    indptr_out[0] = 0;
    int *ptr1, *ptr2, *end1, *end2;
    for (size_t row = 0; row < nrows; row++) {
        
        if (indptr1[row] == indptr1[row+1] ||
            indptr2[row] == indptr2[row+1])
            goto next_row;
        if (ptr_indices1[indptr1[row+1]-1] < ptr_indices2[indptr2[row]] ||
            ptr_indices2[indptr2[row+1]-1] < ptr_indices1[indptr1[row]])
            goto next_row;
        
        ptr1 = ptr_indices1 + indptr1[row];
        ptr2 = ptr_indices2 + indptr2[row];
        end1 = ptr_indices1 + indptr1[row+1];
        end2 = ptr_indices2 + indptr2[row+1];

        while (true)
        {
            if (ptr1 >= end1 || ptr2 >= end2)
                goto next_row;

            else if (*ptr1 == *ptr2) {
                indices_out.push_back(*ptr1);
                values_out.push_back(values1[ptr1 - ptr_indices1] * values2[ptr2 - ptr_indices2]);
                ptr1++;
                ptr2++;
            }

            else if (*ptr1 > *ptr2) {
                ptr2 = std::lower_bound(ptr2, end2, *ptr1);
            }

            else {
                ptr1 = std::lower_bound(ptr1, end1, *ptr2);
            }
        }

        next_row:
        indptr_out[row+1] = indices_out.size();
    }

    VectorConstructorArgs args;
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &indices_out;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.clear(); indices_out.shrink_to_fit();
    args.as_integer = false; args.from_cpp_vec = true; args.num_vec_from = &values_out;
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

template <class RcppVector>
Rcpp::NumericVector multiply_csr_by_dense_elemwise(Rcpp::IntegerVector indptr,
                                                   Rcpp::IntegerVector indices,
                                                   Rcpp::NumericVector values,
                                                   RcppVector dense_mat)
{
    Rcpp::NumericVector values_out(values.size());
    size_t nrows = indptr.size() - 1;
    for (size_t row = 0; row < nrows; row++)
    {
        #pragma omp simd
        for (int el = indptr[row]; el < indptr[row+1]; el++)
        {
            if (std::is_same<RcppVector, Rcpp::LogicalVector>::value)
                values_out[el] = (dense_mat[row + nrows*(size_t)indices[el]] == NA_LOGICAL)?
                                  NA_REAL : values[el] * (bool)dense_mat[row + nrows*(size_t)indices[el]];
            else if (std::is_same<RcppVector, Rcpp::IntegerVector>::value)
                values_out[el] = (dense_mat[row + nrows*(size_t)indices[el]] == NA_INTEGER)?
                                  NA_REAL : values[el] * dense_mat[row + nrows*(size_t)indices[el]];
            else
                values_out[el] = values[el] * dense_mat[row + nrows*(size_t)indices[el]];
        }
    }
    return values_out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csr_by_dense_elemwise_double(Rcpp::IntegerVector indptr,
                                                          Rcpp::IntegerVector indices,
                                                          Rcpp::NumericVector values,
                                                          Rcpp::NumericVector dense_mat)
{
    return multiply_csr_by_dense_elemwise<Rcpp::NumericVector>(indptr, indices, values, dense_mat);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csr_by_dense_elemwise_int(Rcpp::IntegerVector indptr,
                                                       Rcpp::IntegerVector indices,
                                                       Rcpp::NumericVector values,
                                                       Rcpp::IntegerVector dense_mat)
{
    return multiply_csr_by_dense_elemwise<Rcpp::IntegerVector>(indptr, indices, values, dense_mat);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csr_by_dense_elemwise_bool(Rcpp::IntegerVector indptr,
                                                        Rcpp::IntegerVector indices,
                                                        Rcpp::NumericVector values,
                                                        Rcpp::LogicalVector dense_mat)
{
    return multiply_csr_by_dense_elemwise<Rcpp::LogicalVector>(indptr, indices, values, dense_mat);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List add_csr_elemwise(Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
                            Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
                            Rcpp::NumericVector values1, Rcpp::NumericVector values2,
                            const bool substract)
{
    if (indices1.size() == indices2.size() &&
        INTEGER(indptr1) == INTEGER(indptr2) &&
        INTEGER(indices1) == INTEGER(indices2))
    {

        if (substract && REAL(values1) == REAL(values2))
        {
            return Rcpp::List::create(
                Rcpp::_["indptr"] = Rcpp::IntegerVector(indptr1.size()),
                Rcpp::_["indices"] = Rcpp::IntegerVector(),
                Rcpp::_["values"] = Rcpp::NumericVector()
            );
        }

        Rcpp::NumericVector values_out(values1.size());
        if (!substract)
            #pragma omp simd
            for (int el = 0; el < (int)values1.size(); el++)
                values_out[el] = values1[el] + values2[el];
        else
            #pragma omp simd
            for (int el = 0; el < (int)values1.size(); el++)
                values_out[el] = values1[el] - values2[el];
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr1,
            Rcpp::_["indices"] = indices1,
            Rcpp::_["values"] = values_out
        );
    }

    size_t nrows = indptr1.size() - 1;
    Rcpp::List out = Rcpp::List::create(
        Rcpp::_["indptr"] = Rcpp::IntegerVector(nrows+1)
    );

    int *indptr_out = INTEGER(out["indptr"]);
    std::vector<int> indices_out;
    std::vector<double> values_out;
    indices_out.reserve(std::max(indices1.size(), indices2.size()));
    values_out.reserve(std::max(indices1.size(), indices2.size()));

    int *ptr_indices1 = INTEGER(indices1);
    int *ptr_indices2 = INTEGER(indices2);

    indptr_out[0] = 0;
    int *ptr1, *ptr2, *end1, *end2;
    for (size_t row = 0; row < nrows; row++)
    {

        if (indptr1[row] == indptr1[row+1] &&
            indptr2[row] == indptr2[row+1]) {
            
            goto next_row;

        } else if (indptr1[row] == indptr1[row+1]) {
            
            for (int el = indptr2[row]; el < indptr2[row+1]; el++) {
                indices_out.push_back(indices2[el]);
                values_out.push_back(substract? (-values2[el]) : (values2[el]));
            }
            goto next_row;

        } else if (indptr2[row] == indptr2[row+1]) {
            
            for (int el = indptr1[row]; el < indptr1[row+1]; el++) {
                indices_out.push_back(indices1[el]);
                values_out.push_back(values1[el]);
            }
            goto next_row;

        }
        
        ptr1 = ptr_indices1 + indptr1[row];
        ptr2 = ptr_indices2 + indptr2[row];
        end1 = ptr_indices1 + indptr1[row+1];
        end2 = ptr_indices2 + indptr2[row+1];

        while (true)
        {
            if (ptr1 >= end1 || ptr2 >= end2) {
                
                while (ptr1 < end1) {
                    indices_out.push_back(*ptr1);
                    values_out.push_back(values1[ptr1 - ptr_indices1]);
                    ptr1++;
                }
                while (ptr2 < end2) {
                    indices_out.push_back(*ptr2);
                    values_out.push_back(substract? (-values2[ptr2 - ptr_indices2]) : (values2[ptr2 - ptr_indices2]));
                    ptr2++;
                }
                goto next_row;

            }

            else if (*ptr1 == *ptr2) {

                indices_out.push_back(*ptr1);
                values_out.push_back(values1[ptr1 - ptr_indices1] + (substract? (-values2[ptr2 - ptr_indices2]) : values2[ptr2 - ptr_indices2]));
                ptr1++;
                ptr2++;

            }

            else if (*ptr1 > *ptr2) {

                do {
                    indices_out.push_back(*ptr2);
                    values_out.push_back(substract? (-values2[ptr2 - ptr_indices2]) : (values2[ptr2 - ptr_indices2]));
                    ptr2++;
                } while (ptr2 < end2 && *ptr2 < *ptr1);

            }

            else {

                do {
                    indices_out.push_back(*ptr1);
                    values_out.push_back(values1[ptr1 - ptr_indices1]);
                    ptr1++;
                } while (ptr1 < end1 && *ptr1 < *ptr2);
                
            }
        }

        next_row:
        indptr_out[row+1] = indices_out.size();
    }

    VectorConstructorArgs args;
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &indices_out;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.clear(); indices_out.shrink_to_fit();
    args.as_integer = false; args.from_cpp_vec = true; args.num_vec_from = &values_out;
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csr_by_coo_internal
(
    Rcpp::IntegerVector X_csr_indptr_,
    Rcpp::IntegerVector X_csr_indices_,
    Rcpp::NumericVector X_csr_values_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::NumericVector Y_coo_val
)
{
    size_t nnz_y = Y_coo_row.size();
    std::vector<int> out_coo_row;
    std::vector<int> out_coo_col;
    std::vector<double> out_coo_val;
    out_coo_row.reserve(nnz_y);
    out_coo_col.reserve(nnz_y);
    out_coo_val.reserve(nnz_y);
    double val;

    int *restrict X_csr_indptr = INTEGER(X_csr_indptr_);
    int *restrict X_csr_indices = INTEGER(X_csr_indices_);
    double *restrict X_csr_values = REAL(X_csr_values_);

    for (size_t ix = 0; ix < nnz_y; ix++)
    {
        val = extract_single_val_csr(
            X_csr_indptr,
            X_csr_indices,
            X_csr_values,
            Y_coo_row[ix], Y_coo_col[ix],
            false
        );
        if (ISNAN(val) || val != 0)
        {
            out_coo_row.push_back(Y_coo_row[ix]);
            out_coo_val.push_back(Y_coo_col[ix]);
            out_coo_val.push_back(val * Y_coo_val[ix]);
        }
    }


    Rcpp::List out;
    VectorConstructorArgs args;
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &out_coo_row;
    out["row"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    out_coo_row.clear(); out_coo_row.shrink_to_fit();
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &out_coo_col;
    out["col"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    out_coo_col.clear(); out_coo_col.shrink_to_fit();
    args.as_integer = false; args.from_cpp_vec = true; args.num_vec_from = &out_coo_val;
    out["val"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

template <class RcppMatrix, class InputDType>
Rcpp::List multiply_coo_by_dense
(
    RcppMatrix X_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::NumericVector Y_coo_val
)
{
    InputDType *restrict X;
    if (std::is_same<InputDType, float>::value)
        X = (InputDType*)INTEGER(X_);
    else if (std::is_same<RcppMatrix, Rcpp::LogicalMatrix>::value)
        X = (InputDType*)LOGICAL(X_);
    else if (std::is_same<RcppMatrix, Rcpp::IntegerMatrix>::value)
        X = (InputDType*)INTEGER(X_);
    else
        X = (InputDType*)REAL(X_);

    size_t nrows = X_.nrow();
    size_t nnz_y = Y_coo_row.size();
    Rcpp::NumericVector out_coo_val(nnz_y);

    for (size_t ix = 0; ix < nnz_y; ix++)
    {
        if (std::is_same<InputDType, int>::value)
            out_coo_val[ix] = (X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_val[ix]*nrows] == NA_INTEGER)?
                               (NA_REAL) : (Y_coo_val[ix] * X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_val[ix]*nrows]);
        else if (std::is_same<RcppMatrix, Rcpp::LogicalMatrix>::value)
            out_coo_val[ix] = (X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_val[ix]*nrows] == NA_LOGICAL)?
                               (NA_REAL) : (Y_coo_val[ix] * (bool)X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_val[ix]*nrows]);
        else
            out_coo_val[ix] = Y_coo_val[ix] * X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_val[ix]*nrows];
    }

    /* Note: avoid shallow copies as the indices might get sorted in some situations */
    return Rcpp::List::create(
        Rcpp::_["row"] = Rcpp::IntegerVector(Y_coo_row.begin(), Y_coo_row.end()),
        Rcpp::_["col"] = Rcpp::IntegerVector(Y_coo_col.begin(), Y_coo_col.end()),
        Rcpp::_["val"] = out_coo_val
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_coo_by_dense_numeric
(
    Rcpp::NumericMatrix X_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::NumericVector Y_coo_val
)
{
    return multiply_coo_by_dense<Rcpp::NumericMatrix, double>(
        X_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_coo_by_dense_integer
(
    Rcpp::IntegerMatrix X_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::NumericVector Y_coo_val
)
{
    return multiply_coo_by_dense<Rcpp::IntegerMatrix, int>(
        X_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_coo_by_dense_logical
(
    Rcpp::LogicalMatrix X_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::NumericVector Y_coo_val
)
{
    return multiply_coo_by_dense<Rcpp::LogicalMatrix, int>(
        X_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_coo_by_dense_float32
(
    Rcpp::IntegerMatrix X_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::NumericVector Y_coo_val
)
{
    return multiply_coo_by_dense<Rcpp::IntegerMatrix, float>(
        X_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val
    );
}
