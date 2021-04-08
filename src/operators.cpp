#include "MatrixExtra.h"

/*  R's logic for boolean comparisons:
        NA  & TRUE  = NA
        NA  & FALSE = FALSE
        NA  & NA    = NA
        NA  | TRUE  = TRUE
        NA  | FALSE = NA
        NA  |  NA   = NA
*/
#define R_logical_or(x, y) ((x) == NA_LOGICAL)? \
    ( ((y) == NA_LOGICAL)? NA_LOGICAL : ((y)? true : NA_LOGICAL) ) \
        : \
    ( ((y) == NA_LOGICAL)? ((x)? true : NA_LOGICAL) : ((bool)(x) || (bool)(y)) )
#define R_logical_and(x, y) ((x) == NA_LOGICAL)? \
    ( ((y) == NA_LOGICAL)? NA_LOGICAL : ((y)? NA_LOGICAL : false) ) \
        : \
    ( ((y) == NA_LOGICAL)? ((x)? NA_LOGICAL : false) : ((bool)(x) && (bool)(y)) )
#define R_logical_xor(x, y) (((x) == NA_LOGICAL || (y) == NA_LOGICAL)? NA_LOGICAL : ((bool)(x) != (bool)(y)))

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
    InputDType *restrict X = (InputDType*)X_.begin();

    size_t nrows = X_.nrow();
    size_t nnz_y = Y_coo_row.size();
    RcppVector out_coo_val(nnz_y);

    if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
    {
        if (std::is_same<InputDType, int>::value) {
            for (size_t ix = 0; ix < nnz_y; ix++)
                out_coo_val[ix] = (X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_col[ix]*nrows] == NA_INTEGER)?
                                   (NA_REAL) : (Y_coo_val[ix] * X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_col[ix]*nrows]);
        }


        else if (std::is_same<RcppMatrix, Rcpp::LogicalMatrix>::value) {
            for (size_t ix = 0; ix < nnz_y; ix++)
                out_coo_val[ix] = (X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_col[ix]*nrows] == NA_LOGICAL)?
                                   (NA_REAL) : (Y_coo_val[ix] * (bool)X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_col[ix]*nrows]);
        }

        else {
            for (size_t ix = 0; ix < nnz_y; ix++)
                out_coo_val[ix] = Y_coo_val[ix] * X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_col[ix]*nrows];
        }
    }

    else
    {
        for (size_t ix = 0; ix < nnz_y; ix++)
            out_coo_val[ix] = R_logical_and(Y_coo_val[ix], X[(size_t)Y_coo_row[ix] + (size_t)Y_coo_col[ix]*nrows]);
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

template <class RcppMatrix, class InputDType, class DenseType>
Rcpp::List add_NAs_from_dense_after_elemenwise_mult_template
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppMatrix dense_
)
{
    DenseType *restrict dense = (DenseType*)dense_.begin();
    int *restrict indices_begin = indices.begin();
    std::vector<int> ii;
    std::vector<int> jj;
    std::vector<InputDType> xx;

    size_t nrows = dense_.nrow();
    size_t ncols = dense_.ncol();

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (size_t row = 0; row < nrows; row++)
            {
                if (ISNAN(dense[row + col * nrows]))
                {
                    if (indptr[row] == indptr[row+1] ||
                        col < indices_begin[indptr[row]] ||
                        col > indices_begin[indptr[row+1]-1] ||
                        std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], col)
                            >=
                        indices_begin + indptr[row+1])
                    {
                        ii.push_back(row);
                        jj.push_back(col);
                        xx.push_back(NA_REAL);
                    }
                }
            }
        }
    }

    else if (std::is_same<DenseType, float>::value) /* <- float */
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (size_t row = 0; row < nrows; row++)
            {
                if (std::isnan(dense[row + col * nrows]))
                {
                    if (std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], col)
                            ==
                        indices_begin + indptr[row+1])
                    {
                        ii.push_back(row);
                        jj.push_back(col);
                        xx.push_back(NA_REAL);
                    }
                }
            }
        }
    }

    else if (std::is_same<RcppMatrix, Rcpp::IntegerMatrix>::value) /* <- float */
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (size_t row = 0; row < nrows; row++)
            {
                if (dense[row + col * nrows] == NA_INTEGER)
                {
                    if (std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], col)
                            ==
                        indices_begin + indptr[row+1])
                    {
                        ii.push_back(row);
                        jj.push_back(col);
                        xx.push_back(NA_REAL);
                    }
                }
            }
        }
    }

    else
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (size_t row = 0; row < nrows; row++)
            {
                if (dense[row + col * nrows] == NA_LOGICAL)
                {
                    if (std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], col)
                            ==
                        indices_begin + indptr[row+1])
                    {
                        ii.push_back(row);
                        jj.push_back(col);
                        xx.push_back(NA_LOGICAL);
                    }
                }
            }
        }
    }

    Rcpp::List out;
    if (!ii.size())
        return out;

    VectorConstructorArgs args;
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &ii;
    out["ii"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    ii.clear(); ii.shrink_to_fit();
    args.int_vec_from = &jj;
    out["jj"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    jj.clear(); jj.shrink_to_fit();
    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value ||
        std::is_same<DenseType, float>::value
    ) {
        args.as_integer = false; args.num_vec_from = &xx;
    } else {
        args.as_integer = true, args.as_logical = true; args.int_vec_from = &xx;
    }
    out["xx"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List add_NAs_from_dense_after_elemenwise_mult_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericMatrix dense_
)
{
    return add_NAs_from_dense_after_elemenwise_mult_template<Rcpp::NumericMatrix, double, double>(
        indptr,
        indices,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List add_NAs_from_dense_after_elemenwise_mult_integer
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::IntegerMatrix dense_
)
{
    return add_NAs_from_dense_after_elemenwise_mult_template<Rcpp::IntegerMatrix, double, int>(
        indptr,
        indices,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List add_NAs_from_dense_after_elemenwise_mult_float32
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::IntegerMatrix dense_
)
{
    return add_NAs_from_dense_after_elemenwise_mult_template<Rcpp::IntegerMatrix, double, float>(
        indptr,
        indices,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List add_NAs_from_dense_after_elemenwise_mult_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalMatrix dense_
)
{
    return add_NAs_from_dense_after_elemenwise_mult_template<Rcpp::LogicalMatrix, int, int>(
        indptr,
        indices,
        dense_
    );
}

template <class RcppMatrix, class RcppVector, class InputDType>
RcppVector multiply_csc_by_dense_ignore_NAs
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppVector values,
    RcppMatrix dense_
)
{
    size_t ncols = indptr.size() - 1;
    size_t nrows = dense_.nrow();
    size_t nnz = indices.size();
    RcppVector values_out(nnz);
    InputDType *restrict dense = NULL;
    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value) {
        dense = (InputDType*)REAL(dense_);
    } else if (std::is_same<RcppMatrix, Rcpp::LogicalMatrix>::value) {
        dense = (InputDType*)LOGICAL(dense_);
    } else {
        dense = (InputDType*)INTEGER(dense_);
    }

    size_t curr = 0;

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value ||
        std::is_same<InputDType, float>::value)
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (int ix = indptr[col]; ix < indptr[col+1]; ix++)
            {
                values_out[curr++] = values[ix] * dense[(size_t)indices[ix] + col*nrows];
            }
        }
    }

    else if (!std::is_same<RcppVector, Rcpp::LogicalVector>::value)
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (int ix = indptr[col]; ix < indptr[col+1]; ix++)
            {
                values_out[curr++] = (dense[(size_t)indices[ix] + col*nrows] == NA_INTEGER)?
                                      NA_REAL : (values[ix] * dense[(size_t)indices[ix] + col*nrows]);
            }
        }
    }

    else /* <- this is '&' */
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (int ix = indptr[col]; ix < indptr[col+1]; ix++)
            {
                values_out[curr++] = R_logical_and(values[ix], dense[(size_t)indices[ix] + col*nrows]);
            }
        }
    }

    return values_out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csc_by_dense_ignore_NAs_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::NumericMatrix dense_
)
{
    return multiply_csc_by_dense_ignore_NAs<Rcpp::NumericMatrix, Rcpp::NumericVector, double>(
        indptr,
        indices,
        values,
        dense_
    );
}
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csc_by_dense_ignore_NAs_float32
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerMatrix dense_
)
{
    return multiply_csc_by_dense_ignore_NAs<Rcpp::IntegerMatrix, Rcpp::NumericVector, float>(
        indptr,
        indices,
        values,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csc_by_dense_ignore_NAs_integer
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerMatrix dense_
)
{
    return multiply_csc_by_dense_ignore_NAs<Rcpp::IntegerMatrix, Rcpp::NumericVector, int>(
        indptr,
        indices,
        values,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csc_by_dense_ignore_NAs_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::LogicalMatrix dense_
)
{
    return multiply_csc_by_dense_ignore_NAs<Rcpp::LogicalMatrix, Rcpp::NumericVector, int>(
        indptr,
        indices,
        values,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector logicaland_csc_by_dense_ignore_NAs
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    Rcpp::LogicalMatrix dense_
)
{
    return multiply_csc_by_dense_ignore_NAs<Rcpp::LogicalMatrix, Rcpp::LogicalVector, int>(
        indptr,
        indices,
        values,
        dense_
    );
}

template <class RcppVector, class RcppMatrix, class InputDType, class OutputDType>
Rcpp::List multiply_csc_by_dense_keep_NAs_template
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    RcppVector values,
    RcppMatrix dense_
)
{
    size_t ncols = indptr.size() - 1;
    size_t nrows = dense_.nrow();
    InputDType *restrict dense = (InputDType*)dense_.begin();
    int *restrict indices = indices_.begin();

    Rcpp::IntegerVector indptr_out(ncols + 1);
    std::vector<int> indices_out;
    std::vector<OutputDType> values_out;
    indices_out.reserve(indices_.size());
    values_out.reserve(values.size());

    int *ptr, *end;
    int row;

    for (size_t col = 0; col < ncols; col++)
    {
        if (indptr[col] == indptr[col+1])
        {
            if (std::is_same<InputDType, double>::value || std::is_same<InputDType, float>::value)
            {
                for (size_t row = 0; row < nrows; row++)
                {
                    if (ISNAN(dense[row + col*nrows]))
                    {
                        indices_out.push_back(row);
                        values_out.push_back(NA_REAL);
                    }
                }
            }

            else if (std::is_same<InputDType, int>::value)
            {
                for (size_t row = 0; row < nrows; row++)
                {
                    if (dense[row + col*nrows] == NA_INTEGER)
                    {
                        indices_out.push_back(row);
                        values_out.push_back(NA_REAL);
                    }
                }
            }

            else
            {
                for (size_t row = 0; row < nrows; row++)
                {
                    if (dense[row + col*nrows] == NA_LOGICAL)
                    {
                        indices_out.push_back(row);
                        values_out.push_back(NA_LOGICAL);
                    }
                }
            }
            goto next_col;
        }

        ptr = indices + indptr[col];
        end = indices + indptr[col+1];
        row = 0;

        while (true)
        {

            if (ptr >= end) {
                if (std::is_same<InputDType, double>::value || std::is_same<InputDType, float>::value)
                {
                    for (; row < (int)nrows; row++) {
                        if (ISNAN(dense[(size_t)row + col*nrows])) {
                            indices_out.push_back(row);
                            values_out.push_back(NA_REAL);
                        }
                    }
                }

                else if (std::is_same<InputDType, int>::value)
                {
                    for (; row < (int)nrows; row++) {
                        if (dense[(size_t)row + col*nrows] == NA_INTEGER) {
                            indices_out.push_back(row);
                            values_out.push_back(NA_REAL);
                        }
                    }
                }

                else
                {
                    for (; row < (int)nrows; row++) {
                        if (dense[(size_t)row + col*nrows] == NA_LOGICAL) {
                            indices_out.push_back(row);
                            values_out.push_back(NA_LOGICAL);
                        }
                    }
                }
                goto next_col;
            }



            else if (*ptr == row) {
                indices_out.push_back(row);
                if (std::is_same<InputDType, double>::value || std::is_same<InputDType, float>::value)
                    values_out.push_back(values[ptr - indices] * dense[(size_t)row + col*nrows]);
                else if (std::is_same<InputDType, int>::value)
                    values_out.push_back((dense[(size_t)row + col*nrows] == NA_INTEGER)?
                                          NA_REAL : values[ptr - indices] * dense[(size_t)row + col*nrows]);
                else
                    values_out.push_back(R_logical_and(values[ptr - indices], dense[(size_t)row + col*nrows]));
                ptr++;
                row++;
            }


            else if (*ptr < row) {
                ptr = std::lower_bound(ptr, end, row);
            }


            else {
                if (std::is_same<InputDType, double>::value || std::is_same<InputDType, float>::value)
                {
                    for (; row < *ptr; row++) {
                        if (ISNAN(dense[(size_t)row + col*nrows])) {
                            indices_out.push_back(row);
                            values_out.push_back(NA_REAL);
                        }
                    }
                }

                else if (std::is_same<InputDType, int>::value)
                {
                    for (; row < *ptr; row++) {
                        if (dense[(size_t)row + col*nrows] == NA_INTEGER) {
                            indices_out.push_back(row);
                            values_out.push_back(NA_REAL);
                        }
                    }
                }

                else
                {
                    for (; row < *ptr; row++) {
                        if (dense[(size_t)row + col*nrows] == NA_LOGICAL) {
                            indices_out.push_back(row);
                            values_out.push_back(NA_LOGICAL);
                        }
                    }
                }
            }
        }


        next_col:
        indptr_out[col+1] = indices_out.size();
    }


    Rcpp::List out;
    out["indptr"] = indptr_out;
    VectorConstructorArgs args;
    args.as_integer = true; args.from_cpp_vec = true; args.int_vec_from = &indices_out;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.clear(); indices_out.shrink_to_fit();
    if (std::is_same<OutputDType, double>::value) {
        args.as_integer = false; args.num_vec_from = &values_out;
    } else {
        args.as_integer = true; args.as_logical = true; args.int_vec_from = &values_out;
    }
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csc_by_dense_keep_NAs_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::NumericVector values,
    Rcpp::NumericMatrix dense_
)
{
    return multiply_csc_by_dense_keep_NAs_template<
            Rcpp::NumericVector, Rcpp::NumericMatrix, double, double>(
        indptr,
        indices_,
        values,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csc_by_dense_keep_NAs_integer
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::NumericVector values,
    Rcpp::IntegerMatrix dense_
)
{
    return multiply_csc_by_dense_keep_NAs_template<
            Rcpp::NumericVector, Rcpp::IntegerMatrix, int, double>(
        indptr,
        indices_,
        values,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csc_by_dense_keep_NAs_logical
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::NumericVector values,
    Rcpp::LogicalMatrix dense_
)
{
    return multiply_csc_by_dense_keep_NAs_template<
            Rcpp::NumericVector, Rcpp::LogicalMatrix, int, double>(
        indptr,
        indices_,
        values,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csc_by_dense_keep_NAs_float32
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::NumericVector values,
    Rcpp::IntegerMatrix dense_
)
{
    return multiply_csc_by_dense_keep_NAs_template<
            Rcpp::NumericVector, Rcpp::IntegerMatrix, float, double>(
        indptr,
        indices_,
        values,
        dense_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List logicaland_csc_by_dense_keep_NAs
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices_,
    Rcpp::LogicalVector values,
    Rcpp::LogicalMatrix dense_
)
{
    return multiply_csc_by_dense_keep_NAs_template<
            Rcpp::LogicalVector, Rcpp::LogicalMatrix, int, int>(
        indptr,
        indices_,
        values,
        dense_
    );
}
