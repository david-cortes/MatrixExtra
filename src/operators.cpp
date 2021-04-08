#include "MatrixExtra.h"

#define R_logical_or(x, y) ((x == NA_LOGICAL || y == NA_LOGICAL)? NA_LOGICAL : ((bool)x || (bool)y))
#define R_logical_and(x, y) ((x == NA_LOGICAL || y == NA_LOGICAL)? NA_LOGICAL : ((bool)x && (bool)y))
#define R_logical_xor(x, y) ((x == NA_LOGICAL || y == NA_LOGICAL)? NA_LOGICAL : ((bool)x != (bool)y))

/* TODO: change some of the vectors to unique_ptr when there is no dynamic reallocation,
   modify also the slicers. This way it saves memory by not being zero-allocated and
   using lazy allocation. */

/* This function does multiplication and ampersand ("&" operator) */
template <class RcppVector=Rcpp::NumericVector, class InputDType=double>
Rcpp::List multiply_csr_elemwise(Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
                                 Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
                                 RcppVector values1, RcppVector values2)
{
    if (indptr1.size() == indptr2.size() &&
        indices1.size() == indices2.size() &&
        INTEGER(indptr1) == INTEGER(indptr2) &&
        INTEGER(indices1) == INTEGER(indices2)
    ) {
        RcppVector values_out(values1.size());

        if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
            #pragma omp simd
            for (int el = 0; el < (int)values1.size(); el++)
                values_out[el] = values1[el] * values2[el];
        }

        else {
            #pragma omp simd
            for (int el = 0; el < (int)values1.size(); el++)
                values_out[el] = R_logical_and(values1[el], values2[el]);
        }
        
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

    auto nnz_min = std::min(indices1.size(), indices2.size());

    int *indptr_out = INTEGER(out["indptr"]);
    std::unique_ptr<int[]> indices_out(new int[nnz_min]);
    std::unique_ptr<InputDType[]> values_out(new InputDType[nnz_min]);

    int *ptr_indices1 = INTEGER(indices1);
    int *ptr_indices2 = INTEGER(indices2);
    size_t curr = 0;

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
                indices_out[curr] = *ptr1;
                if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
                    values_out[curr] = values1[ptr1 - ptr_indices1] * values2[ptr2 - ptr_indices2];
                else
                    values_out[curr] = R_logical_and(values1[ptr1 - ptr_indices1], values2[ptr2 - ptr_indices2]);
                ptr1++;
                ptr2++;
                curr++;
            }

            else if (*ptr1 > *ptr2) {
                ptr2 = std::lower_bound(ptr2, end2, *ptr1);
            }

            else {
                ptr1 = std::lower_bound(ptr1, end1, *ptr2);
            }
        }

        next_row:
        indptr_out[row+1] = curr;
    }

    VectorConstructorArgs args;
    args.as_integer = true; args.from_pointer = true; args.int_pointer_from = indices_out.get();
    args.cpp_lim_size = true; args.size = curr;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.reset();
    if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
        args.as_integer = false; args.from_pointer = true; args.num_pointer_from = values_out.get();
    } else {
        args.as_integer = true; args.from_pointer = true; args.int_pointer_from = values_out.get();
        args.as_logical = true;
    }
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csr_elemwise
(
    Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
    Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
    Rcpp::NumericVector values1, Rcpp::NumericVector values2
)
{
    return multiply_csr_elemwise<Rcpp::NumericVector, double>(
        indptr1, indptr2,
        indices1, indices2,
        values1, values2
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List logicaland_csr_elemwise
(
    Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
    Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
    Rcpp::LogicalVector values1, Rcpp::LogicalVector values2
)
{
    return multiply_csr_elemwise<Rcpp::LogicalVector, int>(
        indptr1, indptr2,
        indices1, indices2,
        values1, values2
    );
}

template <class RcppVector, class RcppMatrixAsVector>
RcppVector multiply_csr_by_dense_elemwise
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppVector values,
    RcppMatrixAsVector dense_mat
)
{
    RcppVector values_out(values.size());
    size_t nrows = indptr.size() - 1;
    if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
    {
        for (size_t row = 0; row < nrows; row++)
        {
            #pragma omp simd
            for (int el = indptr[row]; el < indptr[row+1]; el++)
            {
                if (std::is_same<RcppMatrixAsVector, Rcpp::LogicalVector>::value)
                    values_out[el] = (dense_mat[row + nrows*(size_t)indices[el]] == NA_LOGICAL)?
                                      NA_REAL : values[el] * (bool)dense_mat[row + nrows*(size_t)indices[el]];
                else if (std::is_same<RcppMatrixAsVector, Rcpp::IntegerVector>::value)
                    values_out[el] = (dense_mat[row + nrows*(size_t)indices[el]] == NA_INTEGER)?
                                      NA_REAL : values[el] * dense_mat[row + nrows*(size_t)indices[el]];
                else
                    values_out[el] = values[el] * dense_mat[row + nrows*(size_t)indices[el]];
            }
        }
    }

    else
    {
        for (size_t row = 0; row < nrows; row++)
        {
            #pragma omp simd
            for (int el = indptr[row]; el < indptr[row+1]; el++)
            {
                values_out[el] = R_logical_and(values[el], dense_mat[row + nrows*(size_t)indices[el]]);
            }
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
    return multiply_csr_by_dense_elemwise(indptr, indices, values, dense_mat);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csr_by_dense_elemwise_float32(Rcpp::IntegerVector indptr,
                                                          Rcpp::IntegerVector indices,
                                                          Rcpp::NumericVector values,
                                                          Rcpp::IntegerVector dense_mat)
{
    return multiply_csr_by_dense_elemwise(indptr, indices, values, (float*)INTEGER(dense_mat));
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csr_by_dense_elemwise_int(Rcpp::IntegerVector indptr,
                                                       Rcpp::IntegerVector indices,
                                                       Rcpp::NumericVector values,
                                                       Rcpp::IntegerVector dense_mat)
{
    return multiply_csr_by_dense_elemwise(indptr, indices, values, dense_mat);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csr_by_dense_elemwise_bool(Rcpp::IntegerVector indptr,
                                                        Rcpp::IntegerVector indices,
                                                        Rcpp::NumericVector values,
                                                        Rcpp::LogicalVector dense_mat)
{
    return multiply_csr_by_dense_elemwise(indptr, indices, values, dense_mat);
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector logicaland_csr_by_dense_cpp
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    Rcpp::LogicalVector dense_mat
)
{
    return multiply_csr_by_dense_elemwise(indptr, indices, values, dense_mat);
}

/* Does "+" and logical "|" */
template <class RcppVector, class InputDType>
Rcpp::List add_csr_elemwise(Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
                            Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
                            RcppVector values1_, RcppVector values2_,
                            const bool substract, const bool xor_op)
{
    if (indices1.size() == indices2.size() &&
        INTEGER(indptr1) == INTEGER(indptr2) &&
        INTEGER(indices1) == INTEGER(indices2))
    {

        if (substract && REAL(values1_) == REAL(values2_))
        {
            return Rcpp::List::create(
                Rcpp::_["indptr"] = Rcpp::IntegerVector(indptr1.size()),
                Rcpp::_["indices"] = Rcpp::IntegerVector(),
                Rcpp::_["values"] = Rcpp::NumericVector()
            );
        }

        RcppVector values_out(values1_.size());
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            if (!substract)
                #pragma omp simd
                for (int el = 0; el < (int)values1_.size(); el++)
                    values_out[el] = values1_[el] + values2_[el];
            else
                #pragma omp simd
                for (int el = 0; el < (int)values1_.size(); el++)
                    values_out[el] = values1_[el] - values2_[el];
        }

        else
        {
            if (!xor_op)
                #pragma omp simd
                for (int el = 0; el < (int)values1_.size(); el++)
                    values_out[el] = R_logical_or(values1_[el], values2_[el]);
            else
                #pragma omp simd
                for (int el = 0; el < (int)values1_.size(); el++)
                    values_out[el] = R_logical_xor(values1_[el], values2_[el]);
        }

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

    size_t nnz_max = (size_t)indices1.size() + (size_t)indices2.size();

    int *indptr_out = INTEGER(out["indptr"]);
    std::unique_ptr<int[]> indices_out(new int[nnz_max]);
    std::unique_ptr<InputDType[]> values_out(new InputDType[nnz_max]);

    int *restrict ptr_indices1 = INTEGER(indices1);
    int *restrict ptr_indices2 = INTEGER(indices2);
    InputDType *restrict values1 = values1_.begin();
    InputDType *restrict values2 = values2_.begin();
    

    indptr_out[0] = 0;
    size_t curr = 0;
    int *ptr1, *ptr2, *end1, *end2;
    for (size_t row = 0; row < nrows; row++)
    {

        if (indptr1[row] == indptr1[row+1] &&
            indptr2[row] == indptr2[row+1]) {
            
            goto next_row;

        } else if (indptr1[row] == indptr1[row+1]) {

            std::copy(ptr_indices2 + indptr2[row], ptr_indices2 + indptr2[row+1], indices_out.get() + curr);
            if (!substract) {
                std::copy(values2 + indptr2[row], values2 + indptr2[row+1], values_out.get() + curr);
                curr += indptr2[row+1] - indptr2[row];
            }
            else {
                for (int el = indptr2[row]; el < indptr2[row+1]; el++)
                    values_out[curr++] = -values2[el];
            }
            
            goto next_row;

        } else if (indptr2[row] == indptr2[row+1]) {
            
            std::copy(ptr_indices1 + indptr1[row], ptr_indices1 + indptr1[row+1], indices_out.get() + curr);
            std::copy(values1 + indptr1[row], values1 + indptr1[row+1], values_out.get() + curr);
            curr += indptr1[row+1] - indptr1[row];
            goto next_row;

        }
        
        ptr1 = ptr_indices1 + indptr1[row];
        ptr2 = ptr_indices2 + indptr2[row];
        end1 = ptr_indices1 + indptr1[row+1];
        end2 = ptr_indices2 + indptr2[row+1];

        while (true)
        {
            if (ptr1 >= end1 || ptr2 >= end2) {
                
                if (ptr1 < end1) {
                    std::copy(ptr1, end1, indices_out.get() + curr);
                    std::copy(values1 + (ptr1 - ptr_indices1), values1 + (end1 - ptr_indices1), values_out.get() + curr);
                    curr += end1 - ptr1;
                }
                if (ptr2 < end2) {
                    std::copy(ptr2, end2, indices_out.get() + curr);
                    if (!substract) {
                        std::copy(values2 + (ptr2 - ptr_indices2), values2 + (end2 - ptr_indices2), values_out.get() + curr);
                        curr += end2 - ptr2;
                    }
                    else {
                        for (int *ix = ptr2; ix < end2; ix++)
                            values_out[curr++] = -values2[ix - ptr_indices2];
                    }
                }
                goto next_row;

            }

            else if (*ptr1 == *ptr2) {

                indices_out[curr] = *ptr1;
                if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
                    values_out[curr++] = values1[ptr1 - ptr_indices1]
                                            + (substract?
                                                (-values2[ptr2 - ptr_indices2])
                                                    :
                                                values2[ptr2 - ptr_indices2]);
                else
                    values_out[curr++] = xor_op?
                                            R_logical_xor(values1[ptr1 - ptr_indices1], values2[ptr2 - ptr_indices2])
                                                :
                                            R_logical_or(values1[ptr1 - ptr_indices1], values2[ptr2 - ptr_indices2]);

                ptr1++;
                ptr2++;

            }

            else if (*ptr1 > *ptr2) {

                do {
                    indices_out[curr] = *ptr2;
                    values_out[curr] = substract? (-values2[ptr2 - ptr_indices2]) : (values2[ptr2 - ptr_indices2]);
                    ptr2++;
                    curr++;
                } while (ptr2 < end2 && *ptr2 < *ptr1);

            }

            else {

                do {
                    indices_out[curr] = *ptr1;
                    values_out[curr] = values1[ptr1 - ptr_indices1];
                    ptr1++;
                    curr++;
                } while (ptr1 < end1 && *ptr1 < *ptr2);
                
            }
        }

        next_row:
        indptr_out[row+1] = curr;
    }

    VectorConstructorArgs args;
    args.as_integer = true; args.from_pointer = true; args.int_pointer_from = indices_out.get();
    args.cpp_lim_size = true; args.size = curr;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.reset();
    if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
        args.as_integer = false; args.from_pointer = true; args.num_pointer_from = values_out.get();
    } else {
        args.as_integer = true; args.from_pointer = true; args.int_pointer_from = values_out.get();
        args.as_logical = true;
    }
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List add_csr_elemwise
(
    Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
    Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
    Rcpp::NumericVector values1, Rcpp::NumericVector values2,
    const bool substract
)
{
    return add_csr_elemwise<Rcpp::NumericVector, double>(
        indptr1, indptr2,
        indices1, indices2,
        values1, values2,
        substract, false
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List logicalor_csr_elemwise
(
    Rcpp::IntegerVector indptr1, Rcpp::IntegerVector indptr2,
    Rcpp::IntegerVector indices1, Rcpp::IntegerVector indices2,
    Rcpp::LogicalVector values1, Rcpp::LogicalVector values2,
    const bool xor_op
)
{
    return add_csr_elemwise<Rcpp::LogicalVector, int>(
        indptr1, indptr2,
        indices1, indices2,
        values1, values2,
        false, xor_op
    );
}

template <class RcppVector, class InputDType>
Rcpp::List multiply_csr_by_coo_elemwise
(
    Rcpp::IntegerVector X_csr_indptr_,
    Rcpp::IntegerVector X_csr_indices_,
    RcppVector X_csr_values_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    RcppVector Y_coo_val,
    int max_row_X, int max_col_X
)
{
    size_t nnz_y = Y_coo_row.size();
    std::unique_ptr<int[]> out_coo_row(new int[nnz_y]);
    std::unique_ptr<int[]> out_coo_col(new int[nnz_y]);
    std::unique_ptr<InputDType[]> out_coo_val(new InputDType[nnz_y]);
    size_t curr = 0;
    InputDType val;

    int *restrict X_csr_indptr = INTEGER(X_csr_indptr_);
    int *restrict X_csr_indices = INTEGER(X_csr_indices_);
    InputDType *restrict X_csr_values = NULL;
    if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
        X_csr_values = (InputDType*)REAL(X_csr_values_);
    } else {
        X_csr_values = (InputDType*)LOGICAL(X_csr_values_);
    }

    if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
    {
        for (size_t ix = 0; ix < nnz_y; ix++)
        {
            if ((ISNAN(Y_coo_val[ix]) || Y_coo_val[ix] != 0) &&
                Y_coo_row[ix] < max_row_X &&
                Y_coo_col[ix] < max_col_X)
            {
                val = extract_single_val_csr(
                    X_csr_indptr,
                    X_csr_indices,
                    X_csr_values,
                    Y_coo_row[ix], Y_coo_col[ix],
                    true
                );

                if (ISNAN(val) || val != 0)
                {
                    out_coo_row[curr] = Y_coo_row[ix];
                    out_coo_col[curr] = Y_coo_col[ix];
                    out_coo_val[curr] = val * Y_coo_val[ix];
                    curr++;
                }
            }
        }
    }

    else
    {
        for (size_t ix = 0; ix < nnz_y; ix++)
        {
            if (Y_coo_val[ix] != 0 &&
                Y_coo_row[ix] < max_row_X &&
                Y_coo_col[ix] < max_col_X)
            {
                val = extract_single_val_csr(
                    X_csr_indptr,
                    X_csr_indices,
                    X_csr_values,
                    Y_coo_row[ix], Y_coo_col[ix],
                    true
                );

                if (val != 0)
                {
                    out_coo_row[curr] = Y_coo_row[ix];
                    out_coo_col[curr] = Y_coo_col[ix];
                    out_coo_val[curr] = R_logical_and(val, Y_coo_val[ix]);
                    curr++;
                }
            }
        }
    }


    Rcpp::List out;
    VectorConstructorArgs args;
    args.as_integer = true; args.from_pointer = true; args.int_pointer_from = out_coo_row.get();
    args.cpp_lim_size = true; args.size = curr;
    out["row"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    out_coo_row.reset();
    args.as_integer = true; args.from_pointer = true; args.int_pointer_from = out_coo_col.get();
    out["col"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    out_coo_col.reset();
    if (std::is_same<RcppVector, Rcpp::NumericVector>::value) {
        args.as_integer = false; args.from_pointer = true; args.num_pointer_from = out_coo_val.get();
    } else {
        args.as_integer = true; args.from_pointer = true; args.int_pointer_from = out_coo_val.get();
        args.as_logical = true;
    }
    out["val"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csr_by_coo_elemwise
(
    Rcpp::IntegerVector X_csr_indptr_,
    Rcpp::IntegerVector X_csr_indices_,
    Rcpp::NumericVector X_csr_values_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::NumericVector Y_coo_val,
    int max_row_X, int max_col_X
)
{
    return multiply_csr_by_coo_elemwise<Rcpp::NumericVector, double>(
        X_csr_indptr_,
        X_csr_indices_,
        X_csr_values_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val,
        max_row_X, max_col_X
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List logicaland_csr_by_coo_elemwise
(
    Rcpp::IntegerVector X_csr_indptr_,
    Rcpp::IntegerVector X_csr_indices_,
    Rcpp::LogicalVector X_csr_values_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::LogicalVector Y_coo_val,
    int max_row_X, int max_col_X
)
{
    return multiply_csr_by_coo_elemwise<Rcpp::LogicalVector, int>(
        X_csr_indptr_,
        X_csr_indices_,
        X_csr_values_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val,
        max_row_X, max_col_X
    );
}

template <class RcppVector, class RcppMatrix, class InputDType>
Rcpp::List multiply_coo_by_dense
(
    RcppMatrix X_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    RcppVector Y_coo_val
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
    RcppVector out_coo_val(nnz_y);

    if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
    {
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
    }

    else
    {
        for (size_t ix = 0; ix < nnz_y; ix++)
            out_coo_val[ix] = R_logical_and(Y_coo_val[ix], X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_val[ix]*nrows]);
    }

    /* Note: avoid shallow copies as the indices might get sorted in some situations */
    /* TODO: revisit this and see if it can be done with shallow copies */
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
    return multiply_coo_by_dense<Rcpp::NumericVector, Rcpp::NumericMatrix, double>(
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
    return multiply_coo_by_dense<Rcpp::NumericVector, Rcpp::IntegerMatrix, int>(
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
    return multiply_coo_by_dense<Rcpp::NumericVector, Rcpp::LogicalMatrix, int>(
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
    return multiply_coo_by_dense<Rcpp::NumericVector, Rcpp::IntegerMatrix, float>(
        X_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List logicaland_coo_by_dense_logical
(
    Rcpp::LogicalMatrix X_,
    Rcpp::IntegerVector Y_coo_row,
    Rcpp::IntegerVector Y_coo_col,
    Rcpp::LogicalVector Y_coo_val
)
{
    return multiply_coo_by_dense<Rcpp::LogicalVector, Rcpp::LogicalMatrix, int>(
        X_,
        Y_coo_row,
        Y_coo_col,
        Y_coo_val
    );
}
