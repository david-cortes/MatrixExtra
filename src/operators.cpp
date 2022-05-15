#include "MatrixExtra.h"

#ifdef __clang__
#   pragma clang diagnostic push
#   pragma clang diagnostic ignored "-Wpass-failed"
#endif

/*  R's logic for boolean comparisons:
        NA  & TRUE  = NA
        NA  & FALSE = FALSE
        NA  & NA    = NA
        NA  | TRUE  = TRUE
        NA  | FALSE = NA
        NA  |  NA   = NA
*/
#ifndef __clang__
#define R_logical_or(x, y) ((x) == NA_LOGICAL)? \
    ( ((y) == NA_LOGICAL)? NA_LOGICAL : ((y)? true : NA_LOGICAL) ) \
        : \
    ( ((y) == NA_LOGICAL)? ((x)? true : NA_LOGICAL) : ((bool)(x) || (bool)(y)) )
#define R_logical_and(x, y) ((x) == NA_LOGICAL)? \
    ( ((y) == NA_LOGICAL)? NA_LOGICAL : ((y)? NA_LOGICAL : false) ) \
        : \
    ( ((y) == NA_LOGICAL)? ((x)? NA_LOGICAL : false) : ((bool)(x) && (bool)(y)) )
#define R_logical_xor(x, y) (((x) == NA_LOGICAL || (y) == NA_LOGICAL)? NA_LOGICAL : ((bool)(x) != (bool)(y)))
#else
/* Attempt at solving issues with CRAN checks, might be a compiler bug with clang. */
int R_logical_or(int x, int y)
{
    if (x == NA_LOGICAL)
    {
        if (y == NA_LOGICAL) {
            return NA_LOGICAL;
        }

        else {
            if (y)
                return 1;
            else
                return NA_LOGICAL;
        }
    }

    else
    {
        if (y == NA_LOGICAL) {
            if (x)
                return 1;
            else
                return NA_LOGICAL;
        }

        else {
            return (bool)x || (bool)y;
        }
    }
}
int R_logical_and(int x, int y)
{
    if (x == NA_LOGICAL)
    {
        if (y == NA_LOGICAL) {
            return NA_LOGICAL;
        }
        else {
            if (y)
                return NA_LOGICAL;
            else
                return 0;
        }
    }

    else
    {
        if (y == NA_LOGICAL) {
            if (x)
                return NA_LOGICAL;
            else
                return 0;
        }

        else {
            return (bool)x && (bool)y;
        }
    }
}
int R_logical_xor(int x, int y)
{
    if (x == NA_LOGICAL || y == NA_LOGICAL)
        return NA_LOGICAL;
    else
        return (bool)x != (bool)y;
}
#endif

/* TODO: some of these operations could benefit from adding 'libdivide' when recycling vectors. */

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
            #ifdef _OPENMP
            #pragma omp simd
            #endif
            for (int el = 0; el < (int)values1.size(); el++)
                values_out[el] = values1[el] * values2[el];
        }

        else {
            #ifdef _OPENMP
            #pragma omp simd
            #endif
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
            #ifdef _OPENMP
            #pragma omp simd
            #endif
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
            #ifdef _OPENMP
            #pragma omp simd
            #endif
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
                #ifdef _OPENMP
                #pragma omp simd
                #endif
                for (int el = 0; el < (int)values1_.size(); el++)
                    values_out[el] = values1_[el] + values2_[el];
            else
                #ifdef _OPENMP
                #pragma omp simd
                #endif
                for (int el = 0; el < (int)values1_.size(); el++)
                    values_out[el] = values1_[el] - values2_[el];
        }

        else
        {
            if (!xor_op)
                #ifdef _OPENMP
                #pragma omp simd
                #endif
                for (int el = 0; el < (int)values1_.size(); el++)
                    values_out[el] = R_logical_or(values1_[el], values2_[el]);
            else
                #ifdef _OPENMP
                #pragma omp simd
                #endif
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

    size_large nnz_max = (size_large)indices1.size() + (size_large)indices2.size();

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
    int *res;
    bool add_el;

    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
    {
        for (size_t col = 0; col < ncols; col++)
        {
            for (size_t row = 0; row < nrows; row++)
            {
                if (ISNAN(dense[row + col * nrows]))
                {
                    add_el = indptr[row] == indptr[row+1] ||
                             (int)col < indices_begin[indptr[row]] ||
                             (int)col > indices_begin[indptr[row+1]-1];
                    if (!add_el) {
                        res = std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], (int)col);
                        add_el = res >= indices_begin + indptr[row+1] || *res != (int)col;
                    }
                    if (add_el)
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
                if (std::isnan((float)dense[row + col * nrows]))
                {
                    add_el = indptr[row] == indptr[row+1] ||
                             (int)col < indices_begin[indptr[row]] ||
                             (int)col > indices_begin[indptr[row+1]-1];
                    if (!add_el) {
                        res = std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], (int)col);
                        add_el = res >= indices_begin + indptr[row+1] || *res != (int)col;
                    }
                    if (add_el)
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
                    add_el = indptr[row] == indptr[row+1] ||
                             (int)col < indices_begin[indptr[row]] ||
                             (int)col > indices_begin[indptr[row+1]-1];
                    if (!add_el) {
                        res = std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], (int)col);
                        add_el = res >= indices_begin + indptr[row+1] || *res != (int)col;
                    }
                    if (add_el)
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
                    add_el = indptr[row] == indptr[row+1] ||
                             (int)col < indices_begin[indptr[row]] ||
                             (int)col > indices_begin[indptr[row+1]-1];
                    if (!add_el) {
                        res = std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], (int)col);
                        add_el = res >= indices_begin + indptr[row+1] || *res != (int)col;
                    }
                    if (add_el)
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

#define recyle_pos(row, col, nrows, len) (((size_large)(row) + (size_large)(col)*(size_large)(nrows)) % (size_large)(len))
#define sign(x) ((x) >= 0)
/* Note: strictly speaking, this is the correct function, but in order to match
   with base R, will use base R's own function, which is inefficient */
// define R_intdiv(x, y) (
//     ISNAN(x)? (x) : (
//         ISNAN(y)? y : (
//             (isinf(x) && isinf(y))? NAN : (
//                 (isinf(y) && (x) && sign(x) != sign(y))? (-1) : floor((x) / (y)) 
//             )
//         )
//     ) 
// )
/* https://github.com/wch/r-source/blob/fe82da3baf849fcd3cc7dbc31c6abc72b57aa083/src/main/arithmetic.c#L217 */
/*  */
#if !defined(HAS_LD) || (HAS_LD == 1)
typedef long double LDOUBLE;
#define c_eps LDBL_EPSILON
#else
typedef double LDOUBLE;
#define c_eps DBL_EPSILON
#endif
static inline double R_intdiv(double x1, double x2)
{
    double q = x1 / x2;
    if (x2 == 0.0 || std::fabs(q) * c_eps > 1 || !R_FINITE(q))
    return q;
    if(std::fabs(q) < 1)
    return (q < 0) ? -1
        : ((x1 < 0 && x2 > 0) ||
           (x1 > 0 && x2 < 0) // differing signs
           ? -1 : 0);
    LDOUBLE tmp = (LDOUBLE)x1 - std::floor(q) * (LDOUBLE)x2;
    return (double) (std::floor(q) + std::floor(tmp/x2));
}

/* Note: strictly speaking, this is the correct function, but in order to match
   with base R, will use base R's own function, which is inefficient */
// define R_modulus(x, y) (
//     isinf(x)? NAN : (
//         isinf(y)? 
//             (isinf(x)? NAN : ( (x == 0 || ISNAN(x) || sign(x) == sign(y))? (x) : (y) )  ) 
//                 : 
//             ( (x) - ((y) * floor((x) / (y))) ) 
//     )
// )
/* https://github.com/wch/r-source/blob/fe82da3baf849fcd3cc7dbc31c6abc72b57aa083/src/main/arithmetic.c#L199 */
/*  */
static inline double R_modulus(double x1, double x2)
{
    if (x2 == 0.0) return R_NaN;
    if(std::fabs(x2) * c_eps > 1 && R_FINITE(x1) && std::fabs(x1) <= std::fabs(x2)) {
    return
        (std::fabs(x1) == std::fabs(x2)) ? 0 :
        ((x1 < 0 && x2 > 0) ||
         (x2 < 0 && x1 > 0))
         ? x1+x2  // differing signs
         : x1   ; // "same" signs (incl. 0)
    }
    double q = x1 / x2;
    LDOUBLE tmp = (LDOUBLE)x1 - std::floor(q) * (LDOUBLE)x2;
    return (double) (tmp - std::floor(tmp/x2) * x2);
}
// define R_pow(x, y) (
//     ((x) == 1)? 1 : (
//         ((y) == 0)? 1 : (
//             ((x) < 0 && isinf(y))? NAN : (
//                 ( (x) < 0 && (((y) < 0 && (y) > (-1)) || (y) != floor(y)) )? NAN : pow(x, y) 
//             ) 
//         )
//     )
// )
/* https://github.com/wch/r-source/blob/fe82da3baf849fcd3cc7dbc31c6abc72b57aa083/src/main/arithmetic.c#L231 */
/*  */
// #if defined(_WIN32) && (SIZE_MAX == UINT64_MAX)
// #   define USE_POWL_IN_R_POW 
// #endif
// static inline double R_pow(double x, double y) /* = x ^ y */
// {
//     /* squaring is the most common of the specially handled cases so
//        check for it first. */
//     if(y == 2.0)
//     return x * x;
//     if(x == 1. || y == 0.)
//     return(1.);
//     if(x == 0.) {
//     if(y > 0.) return(0.);
//     else if(y < 0) return(R_PosInf);
//     else return(y); /* NA or NaN, we assert */
//     }
//     if (R_FINITE(x) && R_FINITE(y)) {
//     /* There was a special case for y == 0.5 here, but
//        gcc 4.3.0 -g -O2 mis-compiled it.  Showed up with
//        100^0.5 as 3.162278, example(pbirthday) failed. */
// #ifdef USE_POWL_IN_R_POW
//     // this is used only on 64-bit Windows (so has powl).
//     return powl(x, y);
// #else
//     return pow(x, y);
// #endif
//     }
//     if (ISNAN(x) || ISNAN(y))
//     return(x + y);
//     if(!R_FINITE(x)) {
//     if(x > 0)       /* Inf ^ y */
//         return (y < 0.)? 0. : R_PosInf;
//     else {          /* (-Inf) ^ y */
//         if(R_FINITE(y) && y == floor(y)) /* (-Inf) ^ n */
//         return (y < 0.) ? 0. : (R_modulus(y, 2.) != 0 ? x  : -x);
//     }
//     }
//     if(!R_FINITE(y)) {
//     if(x >= 0) {
//         if(y > 0)       /* y == +Inf */
//         return (x >= 1) ? R_PosInf : 0.;
//         else        /* y == -Inf */
//         return (x < 1) ? R_PosInf : 0.;
//     }
//     }
//     return R_NaN; // all other cases: (-Inf)^{+-Inf, non-int}; (neg)^{+-Inf}
// }
#define throw_internal_err() Rcpp::stop("Internal error. Please file an issue in GitHub.")
#define max5(a, b, c, d, e) std::max(std::max(std::max(a, b), std::max(c, d)), e)

enum Operation {Multiply, PowerTo, Divide, DivRest, IntDiv};

template <class RcppVector>
RcppVector multiply_csr_by_dvec_no_NAs
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    RcppVector values,
    RcppVector dvec,
    const int ncols,
    const bool multiply,
    const bool powerto,
    const bool divide,
    const bool divrest,
    const bool intdiv,
    const bool X_is_LHS
)
{
    Operation op;
    if (multiply)
        op = Multiply;
    else if (powerto)
        op = PowerTo;
    else if (divide)
        op = Divide;
    else if (divrest)
        op = DivRest;
    else if (intdiv)
        op = IntDiv;
    else
        throw_internal_err();

    RcppVector values_out(values.size());
    const int nrows = indptr.size() - 1;
    const size_t dvec_size = dvec.size();
    double val; int vall;

    /* TODO: some of these can be changed to daxpy */

    if (dvec_size == (size_t)nrows)
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            switch(op)
            {
                case Multiply:
                {
                    for (int row = 0; row < nrows; row++)
                        #ifdef _OPENMP
                        #pragma omp simd
                        #endif
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out[ix] = values[ix] * dvec[row];
                    break;
                }

                case Divide:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = values[ix] / val;
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = val / values[ix];
                        }
                    }
                    break;
                }

                case DivRest:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(values[ix], dvec[row]);
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(dvec[row], values[ix]);
                    }
                    break;
                }

                case IntDiv:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(values[ix], dvec[row]);
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(dvec[row], values[ix]);
                    }
                    break;
                }

                case PowerTo:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(values[ix], dvec[row]);
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(dvec[row], values[ix]);
                    }
                    break;
                }
            }
        }

        else
        {
            for (int row = 0; row < nrows; row++)
                #ifdef _OPENMP
                #pragma omp simd
                #endif
                for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                    values_out[ix] = R_logical_and(values[ix], dvec[row]);
        }
    }

    else if ((size_large)dvec_size >= (size_large)nrows * (size_large)ncols)
    {
        const size_t nrows_ = nrows;
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            switch(op)
            {
                case Multiply:
                {
                    for (size_t row = 0; row < (size_t)nrows; row++)
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out[ix] = values[ix] * dvec[row + (size_t)indices[ix]*nrows_];
                    break;
                }

                case Divide:
                {
                    if (X_is_LHS)
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = values[ix] / dvec[row + (size_t)indices[ix]*nrows_];
                    }

                    else
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = dvec[row + (size_t)indices[ix]*nrows_] / values[ix];
                    }
                    break;
                }

                case DivRest:
                {
                    if (X_is_LHS)
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(values[ix], dvec[row + (size_t)indices[ix]*nrows_]);
                    }

                    else
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(dvec[row + (size_t)indices[ix]*nrows_], values[ix]);
                    }
                    break;
                }

                case IntDiv:
                {
                    if (X_is_LHS)
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(values[ix], dvec[row + (size_t)indices[ix]*nrows_]);
                    }

                    else
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(dvec[row + (size_t)indices[ix]*nrows_], values[ix]);
                    }
                    break;
                }

                case PowerTo:
                {
                    if (X_is_LHS)
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(values[ix], dvec[row + (size_t)indices[ix]*nrows_]);
                    }

                    else
                    {
                        for (size_t row = 0; row < (size_t)nrows; row++)
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(dvec[row + (size_t)indices[ix]*nrows_], values[ix]);
                    }
                    break;
                }
            }
        }

        else
        {
            for (size_t row = 0; row < (size_t)nrows; row++)
                for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                    values_out[ix] = R_logical_and(values[ix], dvec[row + (size_t)indices[ix]*nrows_]);
        }
    }

    else if (dvec_size < (size_t)nrows && ((size_t)nrows % dvec_size) == 0)
    {
        const int dvec_size_ = dvec_size;
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            switch(op)
            {
                case Multiply:
                {
                    for (int row = 0; row < nrows; row++)
                    {
                        val = dvec[row % dvec_size_];
                        #ifdef _OPENMP
                        #pragma omp simd
                        #endif
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out[ix] = values[ix] * val;
                    }
                    break;
                }

                case Divide:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = values[ix] / val;
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = val / values[ix];
                        }
                    }
                    break;
                }

                case DivRest:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(values[ix], val);
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(val, values[ix]);
                        }
                    }
                    break;
                }

                case IntDiv:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(values[ix], val);
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(val, values[ix]);
                        }
                    }
                    break;
                }

                case PowerTo:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(values[ix], val);
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            val = dvec[row % dvec_size_];
                            #ifdef _OPENMP
                            #pragma omp simd
                            #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(val, values[ix]);
                        }
                    }
                    break;
                }
            }
        }

        else
        {
            for (int row = 0; row < nrows; row++)
            {
                vall = dvec[row % dvec_size_];
                #ifdef _OPENMP
                #pragma omp simd
                #endif
                for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                    values_out[ix] = R_logical_and(values[ix], vall);
            }
        }
    }

    else
    {
        if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
        {
            switch(op)
            {
                case Multiply:
                {
                    for (int row = 0; row < nrows; row++)
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out[ix] = values[ix] * dvec[recyle_pos(row, indices[ix], nrows, dvec_size)];
                    }
                    break;
                }

                case Divide:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = values[ix] / dvec[recyle_pos(row, indices[ix], nrows, dvec_size)];
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = dvec[recyle_pos(row, indices[ix], nrows, dvec_size)] / values[ix];
                        }
                    }
                    break;
                }

                case DivRest:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(values[ix], dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]);
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_modulus(dvec[recyle_pos(row, indices[ix], nrows, dvec_size)], values[ix]);
                        }
                    }
                    break;
                }

                case IntDiv:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(values[ix], dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]);
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_intdiv(dvec[recyle_pos(row, indices[ix], nrows, dvec_size)], values[ix]);
                        }
                    }
                    break;
                }

                case PowerTo:
                {
                    if (X_is_LHS)
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(values[ix], dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]);
                        }
                    }

                    else
                    {
                        for (int row = 0; row < nrows; row++)
                        {
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out[ix] = R_pow(dvec[recyle_pos(row, indices[ix], nrows, dvec_size)], values[ix]);
                        }
                    }
                    break;
                }
            }
        }

        else
        {
            for (int row = 0; row < nrows; row++)
            {
                for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                    values_out[ix] = R_logical_and(values[ix], dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]);
            }
        }
    }

    return values_out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_csr_by_dvec_no_NAs_numeric
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::NumericVector dvec,
    const int ncols,
    const bool multiply,
    const bool powerto,
    const bool divide,
    const bool divrest,
    const bool intdiv,
    const bool X_is_LHS
)
{
    return multiply_csr_by_dvec_no_NAs(
        indptr,
        indices,
        values,
        dvec,
        ncols,
        multiply,
        powerto,
        divide,
        divrest,
        intdiv,
        X_is_LHS
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector logicaland_csr_by_dvec_internal
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::LogicalVector values,
    Rcpp::LogicalVector dvec,
    const int ncols
)
{
    return multiply_csr_by_dvec_no_NAs(
        indptr,
        indices,
        values,
        dvec,
        ncols,
        true,
        false,
        false,
        false,
        false,
        true
    );
}

void argsort_buffer_NAs
(
    std::vector<int> &rows_na,
    std::vector<int> &cols_na,
    int *restrict argsorted,
    int *restrict temp
)
{
    if (!rows_na.size())
        return;
    std::iota(argsorted, argsorted + rows_na.size(), (int)0);
    std::sort(argsorted, argsorted + rows_na.size(),
              [&rows_na](const int a, const int b)
              {return rows_na[a] < rows_na[b];});
    for (size_t ix = 0; ix < rows_na.size(); ix++)
        temp[ix] = rows_na[argsorted[ix]];
    std::copy(temp, temp + rows_na.size(), rows_na.begin());
    for (size_t ix = 0; ix < rows_na.size(); ix++)
        temp[ix] = cols_na[argsorted[ix]];
    std::copy(temp, temp + rows_na.size(), cols_na.begin());
}

void add_missing_indices_in_loop
(
    int row, int &curr_row, int nrows, bool &check_sorting,
    std::vector<int>::iterator &curr_pos, std::vector<int>::iterator &next_pos,
    int &n_na, int &n_this,
    std::vector<int> &rows_na, std::vector<int> &cols_na,
    std::vector<int> &indices_out, std::vector<double> &values_out,
    double fill_val
)
{
    if (row == curr_row)
    {
        check_sorting = true;
        next_pos = std::lower_bound(curr_pos, rows_na.end(), row+1);
        n_na = next_pos - curr_pos;
        if (n_na)
        {
            std::copy(cols_na.begin() + (curr_pos - rows_na.begin()),
                      cols_na.begin() + (next_pos - rows_na.begin()),
                      std::back_inserter(indices_out));
            std::fill_n(std::back_inserter(values_out), n_na, fill_val);
            n_this += n_na;
        }

        if (next_pos == rows_na.end())
            curr_row = nrows;
        else
            curr_row = *next_pos;
        curr_pos = next_pos;
    }
}

/* indices need to be sorted for this one */
/* TODO: refactor this one */

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csr_by_dvec_with_NAs
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::NumericVector dvec,
    const int ncols,
    const bool multiply,
    const bool powerto,
    const bool divide,
    const bool divrest,
    const bool intdiv,
    const bool X_is_LHS
)
{
    if ((powerto || divide || divrest) && !X_is_LHS)
        throw_internal_err();
    Operation op;
    if (multiply)
        op = Multiply;
    else if (powerto)
        op = PowerTo;
    else if (divide)
        op = Divide;
    else if (divrest)
        op = DivRest;
    else if (intdiv)
        op = IntDiv;
    else
        throw_internal_err();
    const bool is_div = op == Divide || op == DivRest || op == IntDiv;

    Rcpp::IntegerVector indptr_out;
    std::vector<int> indices_out;
    std::vector<double> values_out;
    const int *restrict indices_begin = indices.begin();

    std::vector<int> rows_na;
    std::vector<int> cols_na;
    std::vector<int> rows_nan;
    std::vector<int> cols_nan;
    std::vector<int> rows_ones;
    std::vector<int> cols_ones;
    std::vector<int> rows_inf;
    std::vector<int> cols_inf;

    const int nrows = indptr.size() - 1;
    const size_t dvec_size = dvec.size();
    bool is_fulldense_case = false;

    double val;
    double *row_val;


    if (dvec_size <= (size_t)nrows && ((size_t)nrows % dvec_size) == 0)
    {


        indptr_out = Rcpp::IntegerVector(indptr.size());
        indices_out.reserve(indices.size());
        values_out.reserve(indices.size());

        const int dvec_size_ = dvec_size;

        switch(op)
        {
            case Multiply:
            {
                for (int row = 0; row < nrows; row++)
                {
                    val = dvec[row % dvec_size_];
                    if (!ISNAN(val) && !std::isinf(val))
                    {
                        std::copy(indices_begin + indptr[row], indices_begin + indptr[row+1],
                                  std::back_inserter(indices_out));
                        #ifdef _OPENMP
                        #pragma omp simd
                        #endif
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(values[ix] * val);
                    }

                    else
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, ISNA(val)? NA_REAL : NAN);

                        if (std::isinf(val)) {
                            row_val = values_out.data() + indptr_out[row];
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                row_val[indices[ix]] = values[ix] * val;
                        }
                    }

                    indptr_out[row+1] = indices_out.size();
                }
                break;
            }

            case Divide:
            {
                for (int row = 0; row < nrows; row++)
                {
                    val = dvec[row % dvec_size_];
                    if (val == 0)
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, NAN);

                        row_val = values_out.data() + indptr_out[row];
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            row_val[indices[ix]] = values[ix] / val;
                    }

                    else if (!ISNAN(val))
                    {
                        std::copy(indices_begin + indptr[row], indices_begin + indptr[row+1],
                                  std::back_inserter(indices_out));
                        // #ifdef _OPENMP
                        // #pragma omp simd
                        // #endif
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(values[ix] / val);
                    }

                    else
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, val);
                    }

                    indptr_out[row+1] = indices_out.size();
                }
                break;
            }

            case DivRest:
            {
                for (int row = 0; row < nrows; row++)
                {
                    val = dvec[row % dvec_size_];
                    if (val == 0)
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, NAN);

                        row_val = values_out.data() + indptr_out[row];
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            row_val[indices[ix]] = R_modulus(values[ix], val);
                    }

                    else if (!ISNAN(val))
                    {
                        std::copy(indices_begin + indptr[row], indices_begin + indptr[row+1],
                                  std::back_inserter(indices_out));
                        // #ifdef _OPENMP
                        // #pragma omp simd
                        // #endif
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_modulus(values[ix], val));
                    }

                    else
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, val);
                    }

                    indptr_out[row+1] = indices_out.size();
                }
                break;
            }

            case IntDiv:
            {
                for (int row = 0; row < nrows; row++)
                {
                    val = dvec[row % dvec_size_];
                    if (val == 0)
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, NAN);

                        row_val = values_out.data() + indptr_out[row];
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            row_val[indices[ix]] = R_intdiv(values[ix], val);
                    }

                    else if (ISNAN(val))
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, val);
                    }

                    else
                    {
                        std::copy(indices_begin + indptr[row], indices_begin + indptr[row+1],
                                  std::back_inserter(indices_out));
                        // #ifdef _OPENMP
                        // #pragma omp simd
                        // #endif
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_intdiv(values[ix], val));
                    }

                    indptr_out[row+1] = indices_out.size();
                }
                break;
            }

            case PowerTo:
            {
                for (int row = 0; row < nrows; row++)
                {
                    val = dvec[row % dvec_size_];
                    if (ISNAN(val))
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, val);

                        row_val = values_out.data() + indptr_out[row];
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            row_val[indices[ix]] = R_pow(values[ix], val);
                    }

                    else if (val <= 0)
                    {
                        for (int col = 0; col < ncols; col++) indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, (val == 0)? 1. : HUGE_VAL);

                        row_val = values_out.data() + indptr_out[row];
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            row_val[indices[ix]] = R_pow(values[ix], val);
                    }

                    else
                    {
                        std::copy(indices_begin + indptr[row], indices_begin + indptr[row+1],
                                  std::back_inserter(indices_out));
                        // #ifdef _OPENMP
                        // #pragma omp simd
                        // #endif
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_pow(values[ix], val));
                    }

                    indptr_out[row+1] = indices_out.size();
                }
                break;
            }
        }


    }

    else if ((size_large)dvec_size >= ((size_large)nrows * (size_large)ncols))
    {


        is_fulldense_case = true;
        const size_t nrows_ = nrows;
        const size_t ncols_ = ncols;
        bool add_el;
        const int *res;


        for (size_t row = 0; row < nrows_; row++)
        {
            for (size_t col = 0; col < ncols_; col++)
            {
                if (ISNAN(dvec[row + col*nrows_]) ||
                    ((is_div || powerto) && dvec[row + col*nrows_] == 0) ||
                    (powerto && dvec[row + col*nrows_] < 0) ||
                    (multiply && std::isinf(dvec[row + col*nrows_])))
                {
                    add_el = indptr[row] == indptr[row+1] ||
                             (int)col < indices_begin[indptr[row]] ||
                             (int)col > indices_begin[indptr[row+1]-1];
                    if (!add_el) {
                        res = std::lower_bound(indices_begin + indptr[row], indices_begin + indptr[row+1], (int)col);
                        add_el = res >= indices_begin + indptr[row+1] || *res != (int)col;
                    }
                    if (add_el)
                    {
                        if ((is_div && dvec[row + col*nrows_] == 0) || ISNA(dvec[row + col*nrows_])) {
                            rows_nan.push_back(row);
                            cols_nan.push_back(col);
                        }

                        else if (powerto && dvec[row + col*nrows_] == 0) {
                            rows_ones.push_back(row);
                            cols_ones.push_back(col);
                        }

                        else if (powerto && dvec[row + col*nrows_] < 0) {
                            rows_inf.push_back(row);
                            cols_inf.push_back(col);
                        }

                        else {
                            rows_na.push_back(row);
                            cols_na.push_back(col);
                        }
                    }
                }
            }
        }
        /* Here will continue in the next condition */


    }

    if (!  ((dvec_size <= (size_t)nrows && ((size_t)nrows % dvec_size) == 0)) )
    {
        /* This is a tricky case. What it will do here is:
            - First determine all the elements which would be NA by recycling
              the vector, store them in COO format.
            - Do the normal loop with the ignore-NAs logic, but at the end of
              each row see if it needs to append elements from the pool of NAs,
              and add them if so, sorting the indices for that row along the way. */

        const size_t nrows_ = nrows;
        const size_large ix_max = (size_large)nrows * (size_large)ncols;
        const size_t n_repeats
            =
        std::ceil((long double)((size_large)nrows * (size_large)ncols) / (long double)dvec_size);

        size_t ix_this, row, col;
        bool add_el;
        const int *res;

        if (!is_fulldense_case)
        {
            for (size_t ix = 0; ix < dvec_size; ix++)
            {
                if (ISNAN(dvec[ix]) ||
                    ((is_div || powerto) && dvec[ix] == 0) ||
                    (powerto && dvec[ix] < 0) ||
                    (multiply && std::isinf(dvec[ix])))
                {
                    for (size_t rep = 0; rep < n_repeats; rep++)
                    {
                        ix_this = ix + rep * dvec_size;
                        if ((size_large)ix_this >= ix_max)
                            break;
                        row = ix_this % nrows_;
                        col = ix_this / nrows_;
                        add_el = indptr[row] == indptr[row+1] ||
                                 (int)col < indices[indptr[row]] ||
                                 (int)col > indices[indptr[row+1]-1];
                        if (!add_el) {
                            res = std::lower_bound(indices_begin + indptr[row],
                                                   indices_begin + indptr[row+1],
                                                   (int)col);
                            add_el = res >= indices_begin + indptr[row+1] || *res != (int)col;
                        }
                        if (add_el)
                        {
                            if ((is_div && dvec[ix] == 0) || ISNA(dvec[ix])) {
                                rows_nan.push_back(row);
                                cols_nan.push_back(col);
                            }

                            else if (powerto && dvec[ix] == 0) {
                                rows_ones.push_back(row);
                                cols_ones.push_back(col);
                            }

                            else if (powerto && dvec[ix] < 0) {
                                rows_inf.push_back(row);
                                cols_inf.push_back(col);
                            }

                            else {
                                rows_na.push_back(row);
                                cols_na.push_back(col);
                            }
                        }
                    }
                }
            }
        }

        if (!rows_na.size() && !rows_nan.size() && !rows_ones.size() && !rows_inf.size())
        {
            return Rcpp::List::create(
                Rcpp::_["indptr"] = indptr,
                Rcpp::_["indices"] = indices,
                Rcpp::_["values"] = multiply_csr_by_dvec_no_NAs(indptr, indices, values, dvec, ncols,
                                                                multiply, powerto, divide, divrest, intdiv, X_is_LHS)
            );
        }


        if (((size_large)rows_na.size() + (size_large)rows_nan.size() +
             (size_large)rows_ones.size() + (size_large)rows_inf.size() +
             (size_large)indices.size())
                >= INT_MAX)
        {
            Rcpp::stop("Error: the resulting matrix would have too many entries for a sparse CSR representation (int overflow).");
        }
        Rcpp::checkUserInterrupt();


        indptr_out = Rcpp::IntegerVector(indptr.size());
        indices_out.reserve(indices.size());
        values_out.reserve(indices.size());
        size_t size_buffer = max5(rows_na.size(), rows_nan.size(), rows_ones.size(), rows_inf.size(), (size_t)ncols);
        std::unique_ptr<int[]> argsorted(new int[size_buffer]);
        std::unique_ptr<int[]> temp(new int[size_times_ratio_dbl(size_buffer)]);
        double *temp_dbl = (double*)temp.get();

        argsort_buffer_NAs(rows_na, cols_na, argsorted.get(), temp.get());
        argsort_buffer_NAs(rows_nan, cols_nan, argsorted.get(), temp.get());
        argsort_buffer_NAs(rows_ones, cols_ones, argsorted.get(), temp.get());
        argsort_buffer_NAs(rows_inf, cols_inf, argsorted.get(), temp.get());

        Rcpp::checkUserInterrupt();

        int curr_row = rows_na.size()? rows_na.front() : nrows;
        auto curr_pos = rows_na.begin();
        auto next_pos = curr_pos;

        int curr_rown = rows_nan.size()? rows_nan.front() : nrows;
        auto curr_posn = rows_nan.begin();
        auto next_posn = curr_pos;

        int curr_row1 = rows_ones.size()? rows_ones.front() : nrows;
        auto curr_pos1 = rows_ones.begin();
        auto next_pos1 = curr_pos1;

        int curr_rowi = rows_inf.size()? rows_inf.front() : nrows;
        auto curr_posi = rows_inf.begin();
        auto next_posi = curr_posi;

        int n_na, n_this;
        int *indices_this;
        double *values_this;
        bool check_sorting;

        for (int row = 0; row < nrows; row++)
        {
            check_sorting = false;
            n_this = 0;

            std::copy(indices_begin + indptr[row], indices_begin + indptr[row+1], std::back_inserter(indices_out));
            if (!is_fulldense_case) {
                switch(op)
                {
                    case Multiply:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(values[ix] * dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]);
                        break;
                    }

                    case Divide:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(values[ix] / dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]);
                        break;
                    }

                    case DivRest:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_modulus(values[ix], dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]));
                        break;
                    }

                    case IntDiv:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_intdiv(values[ix], dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]));
                        break;
                    }

                    case PowerTo:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_pow(values[ix], dvec[recyle_pos(row, indices[ix], nrows, dvec_size)]));
                        break;
                    }
                }
            }

            else {
                switch(op)
                {
                    case Multiply:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(values[ix] * dvec[(size_t)row + (size_t)indices[ix]*nrows_]);
                        break;
                    }

                    case Divide:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(values[ix] / dvec[(size_t)row + (size_t)indices[ix]*nrows_]);
                        break;
                    }

                    case DivRest:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_modulus(values[ix], dvec[(size_t)row + (size_t)indices[ix]*nrows_]));
                        break;
                    }

                    case IntDiv:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_intdiv(values[ix], dvec[(size_t)row + (size_t)indices[ix]*nrows_]));
                        break;
                    }

                    case PowerTo:
                    {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            values_out.push_back(R_pow(values[ix], dvec[(size_t)row + (size_t)indices[ix]*nrows_]));
                        break;
                    }
                }
            }

            add_missing_indices_in_loop(
                row, curr_row, nrows, check_sorting,
                curr_pos, next_pos, n_na, n_this,
                rows_na, cols_na,
                indices_out, values_out,
                NA_REAL
            );

            add_missing_indices_in_loop(
                row, curr_rown, nrows, check_sorting,
                curr_posn, next_posn, n_na, n_this,
                rows_nan, cols_nan,
                indices_out, values_out,
                NAN
            );

            add_missing_indices_in_loop(
                row, curr_row1, nrows, check_sorting,
                curr_pos1, next_pos1, n_na, n_this,
                rows_ones, cols_ones,
                indices_out, values_out,
                1.
            );

            add_missing_indices_in_loop(
                row, curr_rowi, nrows, check_sorting,
                curr_posi, next_posi, n_na, n_this,
                rows_inf, cols_inf,
                indices_out, values_out,
                HUGE_VAL
            );

            if (check_sorting)
            {
                indices_this = indices_out.data() + indptr_out[row];
                values_this = values_out.data() + indptr_out[row];
                n_this += indptr[row+1] - indptr[row];

                if (!check_is_sorted(indices_this, n_this))
                {
                    std::iota(argsorted.get(), argsorted.get() + n_this, (size_t)0);
                    std::sort(argsorted.get(), argsorted.get() + n_this,
                              [&indices_this](const int a, const int b)
                              {return indices_this[a] < indices_this[b];});
                    for (int ix = 0; ix < n_this; ix++)
                        temp[ix] = indices_this[argsorted[ix]];
                    std::copy(temp.get(), temp.get() + n_this, indices_this);
                    for (int ix = 0; ix < n_this; ix++)
                        temp_dbl[ix] = values_this[argsorted[ix]];
                    std::copy(temp_dbl, temp_dbl + n_this, values_this);
                }
            }

            indptr_out[row+1] = indices_out.size();
        }
    }

    Rcpp::List out;
    out["indptr"] = indptr_out;
    VectorConstructorArgs args;
    args.from_cpp_vec = true; args.as_integer = true; args.int_vec_from = &indices_out;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.clear(), indices_out.shrink_to_fit();
    args.as_integer = false; args.num_vec_from = &values_out;
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}


template <class RcppVector>
RcppVector multiply_coo_by_dense_ignore_NAs_template
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    RcppVector xx,
    RcppVector dvec,
    const int nrows,
    const int ncols,
    const bool multiply,
    const bool powerto,
    const bool divide,
    const bool divrest,
    const bool intdiv,
    const bool X_is_LHS
)
{
    Operation op;
    if (multiply)
        op = Multiply;
    else if (powerto)
        op = PowerTo;
    else if (divide)
        op = Divide;
    else if (divrest)
        op = DivRest;
    else if (intdiv)
        op = IntDiv;
    else
        throw_internal_err();

    RcppVector xx_out(xx.size());
    const size_t dvec_size = dvec.size();
    const int dvec_size_ = dvec_size;
    const size_t nnz = ii.size();
    if ((size_t)xx_out.size() != nnz)
        Rcpp::stop("Unexpected error.");
    const size_t nrows_ = nrows;

    if (dvec_size == (size_t)nrows)
    {
        switch(op)
        {
            case Multiply:
            {
                if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] * dvec[ii[ix]];
                }

                else
                {
                    // #ifdef _OPENMP
                    // #pragma omp simd
                    // #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_logical_and(xx[ix], dvec[ii[ix]]);
                }
                break;
            }

            case Divide:
            {
                if (X_is_LHS)
                {
                    // #ifdef _OPENMP
                    // #pragma omp simd
                    // #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] / dvec[ii[ix]];
                }

                else
                {
                    // #ifdef _OPENMP
                    // #pragma omp simd
                    // #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = dvec[ii[ix]] / xx[ix];
                }
                break;
            }

            case IntDiv:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(xx[ix], dvec[ii[ix]]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(dvec[ii[ix]], xx[ix]);
                }
                break;
            }

            case DivRest:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(xx[ix], dvec[ii[ix]]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(dvec[ii[ix]], xx[ix]);
                }
                break;
            }

            case PowerTo:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(xx[ix], dvec[ii[ix]]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(dvec[ii[ix]], xx[ix]);
                }
                break;
            }
        }
    }

    else if ((size_large)dvec_size_ == (size_large)nrows * (size_large)ncols)
    {
        switch(op)
        {
            case Multiply:
            {
                if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] * dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_];
                }

                else
                {
                    // #ifdef _OPENMP
                    // #pragma omp simd
                    // #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_logical_and(xx[ix], dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_]);
                }
                break;
            }

            case Divide:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] / dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_];
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_] / xx[ix];
                }
                break;
            }

            case IntDiv:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(xx[ix], dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_], xx[ix]);
                }
                break;
            }

            case DivRest:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(xx[ix], dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_], xx[ix]);
                }
                break;
            }

            case PowerTo:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(xx[ix], dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(dvec[(size_t)ii[ix] + (size_t)jj[ix]*nrows_], xx[ix]);
                }
                break;
            }
        }
    }

    else if (dvec_size < (size_t)nrows && ((size_t)nrows % dvec_size) == 0)
    {
        switch(op)
        {
            case Multiply:
            {
                if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] * dvec[ii[ix] % dvec_size_];
                }

                else
                {
                    // #ifdef _OPENMP
                    // #pragma omp simd
                    // #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_logical_and(xx[ix], dvec[ii[ix] % dvec_size_]);
                }
                break;
            }

            case Divide:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] / dvec[ii[ix] % dvec_size_];
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = dvec[ii[ix] % dvec_size_] / xx[ix];
                }
                break;
            }

            case IntDiv:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(xx[ix], dvec[ii[ix] % dvec_size_]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(dvec[ii[ix] % dvec_size_], xx[ix]);
                }
                break;
            }

            case DivRest:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(xx[ix], dvec[ii[ix] % dvec_size_]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(dvec[ii[ix] % dvec_size_], xx[ix]);
                }
                break;
            }

            case PowerTo:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(xx[ix], dvec[ii[ix] % dvec_size_]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(dvec[ii[ix] % dvec_size_], xx[ix]);
                }
                break;
            }
        }
    }

    else
    {
        switch(op)
        {
            case Multiply:
            {
                if (std::is_same<RcppVector, Rcpp::NumericVector>::value)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] * dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)];
                }

                else
                {
                    // #ifdef _OPENMP
                    // #pragma omp simd
                    // #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_logical_and(xx[ix], dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)]);
                }
                break;
            }

            case Divide:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = xx[ix] / dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)];
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)] / xx[ix];
                }
                break;
            }

            case IntDiv:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(xx[ix], dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_intdiv(dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)], xx[ix]);
                }
                break;
            }

            case DivRest:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(xx[ix], dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_modulus(dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)], xx[ix]);
                }
                break;
            }

            case PowerTo:
            {
                if (X_is_LHS)
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(xx[ix], dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)]);
                }

                else
                {
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (size_t ix = 0; ix < nnz; ix++)
                        xx_out[ix] = R_pow(dvec[recyle_pos(ii[ix], jj[ix], nrows_, dvec_size)], xx[ix]);
                }
                break;
            }
        }
    }


    return xx_out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multiply_coo_by_dense_ignore_NAs_numeric
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::NumericVector xx,
    Rcpp::NumericVector dvec,
    const int nrows,
    const int ncols,
    const bool multiply,
    const bool powerto,
    const bool divide,
    const bool divrest,
    const bool intdiv,
    const bool X_is_LHS
)
{
    return multiply_coo_by_dense_ignore_NAs_template(
        ii,
        jj,
        xx,
        dvec,
        nrows,
        ncols,
        multiply,
        powerto,
        divide,
        divrest,
        intdiv,
        X_is_LHS
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector multiply_coo_by_dense_ignore_NAs_logical
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::LogicalVector xx,
    Rcpp::LogicalVector dvec,
    const int nrows,
    const int ncols
)
{
    return multiply_coo_by_dense_ignore_NAs_template(
        ii,
        jj,
        xx,
        dvec,
        nrows,
        ncols,
        true,
        false,
        false,
        false,
        false,
        true
    );
}


/* For sparse vectors, it should only attempt to do it when the input length is
   a factor of, or equal to, the row size, otherwise it's more efficient in CSC
   and can leave it to 'Matrix' instead */
// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csr_by_svec_no_NAs
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector ii_base1,
    Rcpp::NumericVector xx,
    const int length
)
{
    int nrows = indptr.size() - 1;
    Rcpp::IntegerVector indptr_out(indptr.size());

    if (!ii_base1.size()) {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr_out,
            Rcpp::_["indices"] = Rcpp::IntegerVector(),
            Rcpp::_["values"] = Rcpp::NumericVector()
        );
    }

    std::unique_ptr<int[]> indices_out(new int[indices.size()]);
    std::unique_ptr<double[]> values_out(new double[indices.size()]);
    const bool is_binary = xx.size() == 0;
    const int n_recycles = nrows / length;
    const int nnz_v = ii_base1.size();

    int curr = 0;
    double val;
    int row;

    for (int recycle = 0; recycle < n_recycles; recycle++)
    {
        for (int el = 0; el < nnz_v; el++)
        {
            row = ii_base1[el] - 1 + (recycle * length);
            std::copy(indices.begin() + indptr[row], indices.begin() + indptr[row+1], indices_out.get() + curr);
            indptr_out[row+1] = indptr[row+1] - indptr[row];

            if (!is_binary)
            {
                val = xx[el];
                #ifdef _OPENMP
                #pragma omp simd
                #endif
                for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                    values_out[curr++] = values[ix] * val;
            }

            else
            {
                std::copy(values.begin() + indptr[row], values.begin() + indptr[row+1], values_out.get() + curr);
                curr += indptr_out[row+1];
            }
        }
    }

    for (int row = 0; row < nrows; row++)
        indptr_out[row+1] += indptr_out[row];


    Rcpp::List out;
    out["indptr"] = indptr_out;
    VectorConstructorArgs args;
    args.from_pointer = true; args.as_integer = true; args.int_pointer_from = indices_out.get();
    args.cpp_lim_size = true; args.size = curr;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.reset();
    args.as_integer = false; args.num_pointer_from = values_out.get();
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_csr_by_svec_keep_NAs
(
    Rcpp::IntegerVector indptr,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector ii_base1,
    Rcpp::NumericVector xx,
    const int ncols,
    const int length
)
{
    int nrows = indptr.size() - 1;
    Rcpp::IntegerVector indptr_out(indptr.size());
    const bool X_has_NAs_or_inf = contains_any_nas_or_inf(values);

    if (!ii_base1.size() && !X_has_NAs_or_inf) {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr_out,
            Rcpp::_["indices"] = Rcpp::IntegerVector(),
            Rcpp::_["values"] = Rcpp::NumericVector()
        );
    }

    std::vector<int> indices_out;
    std::vector<double> values_out;
    indices_out.reserve(indices.size());
    values_out.reserve(indices.size());
    const bool is_binary = xx.size() == 0;
    const int n_recycles = nrows / length;
    const int nnz_v = ii_base1.size();

    if (nnz_v > length || length > nrows || (nrows % length) != 0)
        throw_internal_err();

    double val;
    int row;
    int row_end;
    int offset;
    std::vector<double>::iterator this_row;
    const auto ii_begin = ii_base1.begin();
    auto curr_i = ii_base1.begin();
    const auto end_i = ii_base1.end();

    for (int recycle = 0; recycle < n_recycles; recycle++)
    {
        if (X_has_NAs_or_inf)
        {
            curr_i = ii_base1.begin();
            row = recycle * length;
            offset = row;
            row_end = offset + length;

            while (true)
            {
                if (curr_i >= end_i || row >= row_end) {

                    for (; row < row_end; row++) {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++) {
                            if (ISNAN(values[ix]) || std::isinf(values[ix])) {
                                values_out.push_back(std::isinf(values[ix])? NAN : values[ix]);
                                indices_out.push_back(indices[ix]);
                                indptr_out[row+1] = indptr_out[row+1] + 1;
                            }
                        }
                    }

                    break;

                }

                else if (row == *curr_i - 1 + offset) {

                    val = is_binary? 1 : xx[curr_i - ii_begin];

                    if (!is_binary && (ISNAN(val) || std::isinf(val)))
                    {
                        indptr_out[row+1] = ncols;
                        for (int col = 0; col < ncols; col++)
                            indices_out.push_back(col);
                        std::fill_n(std::back_inserter(values_out), ncols, ISNAN(val)? val : NAN);
                        if (std::isinf(val)) {
                            this_row = values_out.end() - ncols;
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                this_row[indices[ix]] = val * values[ix];
                        }
                    }

                    else
                    {
                        indptr_out[row+1] = indptr[row+1] - indptr[row];
                        std::copy(indices.begin() + indptr[row],
                                  indices.begin() + indptr[row+1],
                                  std::back_inserter(indices_out));

                        if (!is_binary) {
                            // #ifdef _OPENMP
                            // #pragma omp simd
                            // #endif
                            for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                                values_out.push_back(values[ix] * val);
                        }

                        else {
                            std::copy(values.begin() + indptr[row],
                                      values.begin() + indptr[row+1],
                                      std::back_inserter(values_out));
                        }
                    }

                    row++;
                    curr_i++;

                }

                else if (row > *curr_i - 1 + offset) {

                    curr_i = std::lower_bound(curr_i, end_i, row + 1 - offset);

                }

                else {

                    for (; row < *curr_i - 1 + offset; row++) {
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++) {
                            if (ISNAN(values[ix]) || std::isinf(values[ix])) {
                                values_out.push_back(std::isinf(values[ix])? NAN : values[ix]);
                                indices_out.push_back(indices[ix]);
                                indptr_out[row+1] = indptr_out[row+1] + 1;
                            }
                        }
                    }

                }
            }
        }

        else
        {
            for (int el = 0; el < nnz_v; el++)
            {
                row = ii_base1[el] - 1 + (recycle * length);
                
                if (!is_binary && (ISNAN(xx[el]) || std::isinf(xx[el])))
                {
                    val = xx[el];
                    indptr_out[row+1] = ncols;
                    for (int col = 0; col < ncols; col++)
                        indices_out.push_back(col);
                    std::fill_n(std::back_inserter(values_out), ncols, ISNAN(val)? val : NAN);
                    if (std::isinf(val)) {
                        this_row = values_out.end() - ncols;
                        for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                            this_row[indices[ix]] = val * values[ix];
                    }

                    continue;
                }

                indptr_out[row+1] = indptr[row+1] - indptr[row];
                std::copy(indices.begin() + indptr[row],
                          indices.begin() + indptr[row+1],
                          std::back_inserter(indices_out));

                if (!is_binary) {
                    val = xx[el];
                    // #ifdef _OPENMP
                    // #pragma omp simd
                    // #endif
                    for (int ix = indptr[row]; ix < indptr[row+1]; ix++)
                        values_out.push_back(values[ix] * val);
                }

                else {
                    std::copy(values.begin() + indptr[row],
                              values.begin() + indptr[row+1],
                              std::back_inserter(values_out));
                }
            }
        }
    }


    for (int row = 0; row < nrows; row++)
        indptr_out[row+1] += indptr_out[row];


    Rcpp::List out;
    out["indptr"] = indptr_out;
    VectorConstructorArgs args;
    args.from_cpp_vec = true; args.as_integer = true; args.int_vec_from = &indices_out;
    args.cpp_lim_size = false;
    out["indices"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    indices_out.clear(); indices_out.shrink_to_fit();
    args.as_integer = false; args.num_vec_from = &values_out;
    out["values"] = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    return out;
}

template <class RcppMatrix, class InputDType>
Rcpp::List multiply_elemwise_dense_by_svec_template
(
    RcppMatrix X_,
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const int length,
    const bool keep_NAs
)
{
    const size_t nnz = ii.size();
    const int nrows = X_.nrow();
    const int ncols = X_.ncol();
    const auto* restrict X = (InputDType*)X_.begin();
    const size_t nrows_ = nrows;
    const size_t ncols_ = ncols;
    size_t row, col;

    const bool simple_mult = std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value ||
                             std::is_same<InputDType, float>::value;

    if ((size_t)length == (size_t)nrows * (size_t)ncols)
    {
        Rcpp::NumericMatrix out_(nrows, ncols);
        auto out = out_.begin();

        if (keep_NAs)
        {
            size_t tot = (size_t)nrows * (size_t)ncols;
            
            if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value) {
                for (size_t ix = 0; ix < tot; ix++)
                    out[ix] = ISNAN(X[ix])? X[ix] : (std::isinf((double)X[ix])? NAN : 0);
            }

            else if (std::is_same<InputDType, float>::value) {
                for (size_t ix = 0; ix < tot; ix++)
                    out[ix] = (std::isnan((float)X[ix]) || std::isinf((float)X[ix])? NAN : 0);
            }

            else {
                for (size_t ix = 0; ix < tot; ix++)
                    out[ix] = (X[ix] == NA_INTEGER)? NA_REAL : 0;
            }
        }

        if (simple_mult)
        {
            for (size_t ix = 0; ix < nnz; ix++) {
                row = (ii[ix]-1) % nrows;
                col = (ii[ix]-1) / nrows;
                out[row + col*nrows_] = X[row + col*nrows_] * xx[ix];
            }
        }

        else
        {
            for (size_t ix = 0; ix < nnz; ix++) {
                row = (ii[ix]-1) % nrows;
                col = (ii[ix]-1) / nrows;
                out[row + col*nrows_] = (X[row + col*nrows_] == NA_INTEGER)? NAN : (X[row + col*nrows_] * xx[ix]);
            }
        }
        return Rcpp::List::create(Rcpp::_["X_dense"] = out_);
    }

    else if (length == nrows)
    {
        Rcpp::IntegerVector indptr(nrows+1);
        Rcpp::IntegerVector indices;
        Rcpp::NumericVector values;
        size_t curr = 0;

        if (!keep_NAs)
        {
            indices = Rcpp::IntegerVector(nnz * (size_t)ncols);
            values = Rcpp::NumericVector(nnz * (size_t)ncols);

            if (simple_mult)
            {
                for (size_t ix = 0; ix < nnz; ix++)
                {
                    row = ii[ix] - 1;
                    indptr[row+1] = ncols;
                    std::iota(indices.begin() + curr, indices.begin() + curr + (size_t)ncols, 0);
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int col = 0; col < ncols; col++)
                        values[curr++] = X[row + col*nrows_] * xx[ix];
                }
            }

            else
            {
                for (size_t ix = 0; ix < nnz; ix++)
                {
                    row = ii[ix] - 1;
                    indptr[row+1] = ncols;
                    std::iota(indices.begin() + curr, indices.begin() + curr + (size_t)ncols, 0);
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int col = 0; col < ncols; col++)
                        values[curr++] = (X[row + col*nrows_] == NA_INTEGER)? NA_REAL : (X[row + col*nrows_] * xx[ix]);
                }
            }

            for (int row = 0; row < nrows; row++)
                indptr[row+1] += indptr[row];
        }

        else
        {
            std::vector<int> indices_;
            std::vector<double> values_;
            indices_.reserve(nnz * (size_t)ncols);
            values_.reserve(nnz * (size_t)ncols);

            const auto ii_begin = ii.begin();
            auto curr_i = ii_begin;
            const auto end_i = ii.end();
            int row = 0;

            while (true)
            {
                if (curr_i >= end_i || row >= nrows) {

                    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
                    {
                        for (; row < nrows; row++) {
                            for (size_t col = 0; col < ncols_; col++) {
                                if (ISNAN(X[(size_t)row + col*nrows_]) ||
                                    std::isinf((double)X[(size_t)row + col*nrows_])
                                ) {
                                    indices_.push_back(col);
                                    values_.push_back(
                                        std::isinf((double)X[(size_t)row + col*nrows_])?
                                            NAN : X[(size_t)row + col*nrows_]
                                    );
                                }
                            }
                            indptr[row+1] = indices_.size();
                        }
                    }

                    else if (std::is_same<InputDType, float>::value)
                    {
                        for (; row < nrows; row++) {
                            for (size_t col = 0; col < ncols_; col++) {
                                if (std::isnan((float)X[(size_t)row + col*nrows_]) ||
                                    std::isinf((float)X[(size_t)row + col*nrows_])
                                ) {
                                    indices_.push_back(col);
                                    values_.push_back(
                                        NAN
                                    );
                                }
                            }
                            indptr[row+1] = indices_.size();
                        }
                    }

                    else
                    {
                        for (; row < nrows; row++) {
                            for (size_t col = 0; col < ncols_; col++) {
                                if (X[(size_t)row + col*nrows_] == NA_INTEGER
                                ) {
                                    indices_.push_back(col);
                                    values_.push_back(
                                        NA_REAL
                                    );
                                }
                            }
                            indptr[row+1] = indices_.size();
                        }
                    }

                    break;

                }

                else if (row == *curr_i-1) {

                    if (simple_mult)
                    {
                        for (size_t col = 0; col < ncols_; col++) {
                            indices_.push_back(col);
                            values_.push_back(X[(size_t)row + col*nrows_] * xx[curr_i - ii_begin]);
                        }
                    }

                    else
                    {
                        for (size_t col = 0; col < ncols_; col++) {
                            indices_.push_back(col);
                            values_.push_back(
                                (X[(size_t)row + col*nrows_] == NA_INTEGER)? NAN :
                                (X[(size_t)row + col*nrows_] * xx[curr_i - ii_begin])
                            );
                        }
                    }
                    indptr[row+1] = indices_.size();
                    row++;
                    curr_i++;

                }

                else if (row > *curr_i-1) {

                    curr_i = std::lower_bound(curr_i, end_i, row+1);

                }

                else {

                    if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
                    {
                        for (; row < *curr_i-1; row++) {
                            for (size_t col = 0; col < ncols_; col++) {
                                if (ISNAN(X[(size_t)row + col*nrows_]) ||
                                    std::isinf((double)X[(size_t)row + col*nrows_])
                                ) {
                                    indices_.push_back(col);
                                    values_.push_back(
                                        std::isinf((double)X[(size_t)row + col*nrows_])?
                                            NAN : X[(size_t)row + col*nrows_]
                                    );
                                }
                            }
                            indptr[row+1] = indices_.size();
                        }
                    }

                    else if (std::is_same<InputDType, float>::value)
                    {
                        for (; row < *curr_i-1; row++) {
                            for (size_t col = 0; col < ncols_; col++) {
                                if (std::isnan((float)X[(size_t)row + col*nrows_]) ||
                                    std::isinf((float)X[(size_t)row + col*nrows_])
                                ) {
                                    indices_.push_back(col);
                                    values_.push_back(
                                        NAN
                                    );
                                }
                            }
                            indptr[row+1] = indices_.size();
                        }
                    }

                    else
                    {
                        for (; row < *curr_i-1; row++) {
                            for (size_t col = 0; col < ncols_; col++) {
                                if (X[(size_t)row + col*nrows_] == NA_INTEGER) {
                                    indices_.push_back(col);
                                    values_.push_back(
                                        NA_REAL
                                    );
                                }
                            }
                            indptr[row+1] = indices_.size();
                        }
                    }

                }
            }

            VectorConstructorArgs args;
            args.from_cpp_vec = true; args.as_integer = true; args.int_vec_from = &indices_;
            indices = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
            indices_.clear(); indices_.shrink_to_fit();
            args.as_integer = false; args.num_vec_from = &values_;
            values = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
        }

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    else if (length < nrows && (nrows % length) == 0)
    {
        size_t curr = 0;
        int n_recycles = nrows / length;
        const int one = 1;
        double *restrict xx_ = REAL(xx);

        Rcpp::IntegerVector indptr(nrows+1);
        Rcpp::IntegerVector indices;
        Rcpp::NumericVector values;

        if (!keep_NAs)
        {
            indices = Rcpp::IntegerVector((size_t)n_recycles * nnz * ncols_);
            values = Rcpp::NumericVector((size_t)n_recycles * nnz * ncols_);
            double *restrict values_ = REAL(values);
            double *restrict X__ = (double*)X;
            double val;

            for (int recycle = 0; recycle < n_recycles; recycle++)
            {
                if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
                {
                    for (size_t ix = 0; ix < nnz; ix++)
                    {
                        row = ii[ix] - 1 + recycle*length;
                        indptr[row+1] = ncols;
                        std::iota(indices.begin() + curr, indices.begin() + curr + ncols_, 0);
                        daxpy_(&ncols, xx_ + ix, X__ + row, &nrows, values_ + curr, &one);
                        curr += ncols;
                    }
                }

                else if (std::is_same<InputDType, float>::value)
                {
                    for (size_t ix = 0; ix < nnz; ix++)
                    {
                        row = ii[ix] - 1 + recycle*length;
                        val = xx_[ix];
                        indptr[row+1] = ncols;
                        std::iota(indices.begin() + curr, indices.begin() + curr + ncols_, 0);
                        for (size_t col = 0; col < ncols_; col++)
                            values[curr++] = X[(size_t)row + col*nrows_] * val;
                    }
                }

                else
                {
                    for (size_t ix = 0; ix < nnz; ix++)
                    {
                        row = ii[ix] - 1 + recycle*length;
                        val = xx_[ix];
                        indptr[row+1] = ncols;
                        std::iota(indices.begin() + curr, indices.begin() + curr + ncols_, 0);
                        for (size_t col = 0; col < ncols_; col++)
                            values[curr++] = (X[(size_t)row + col*nrows_] == NA_INTEGER)? NAN : (X[(size_t)row + col*nrows_] * val);
                    }
                }
            }

            for (int row = 0; row < nrows; row++)
                indptr[row+1] += indptr[row];
        }

        else
        {
            std::vector<int> indices_;
            std::vector<double> values_;
            indices_.reserve((size_t)n_recycles * nnz * ncols_);
            values_.reserve((size_t)n_recycles * nnz * ncols_);

            int row = 0;
            int row_end;
            int offset;

            const auto ii_begin = ii.begin();
            auto curr_i = ii_begin;
            const auto end_i = ii.end();

            for (int recycle = 0; recycle < n_recycles; recycle++)
            {
                offset = recycle * length;
                row = offset;
                row_end = row + length;
                
                while (true)
                {
                    if (curr_i >= end_i || row >= row_end) {

                        if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
                        {
                            for (; row < row_end; row++) {
                                for (size_t col = 0; col < ncols_; col++) {
                                    if (ISNAN(X[(size_t)row + col*nrows_]) ||
                                        std::isinf((double)X[(size_t)row + col*nrows_])
                                    ) {
                                        indices_.push_back(col);
                                        values_.push_back(
                                            std::isinf((double)X[(size_t)row + col*nrows_])?
                                                NAN : X[(size_t)row + col*nrows_]
                                        );
                                    }
                                }
                                indptr[row+1] = indices_.size();
                            }
                        }

                        else if (std::is_same<InputDType, float>::value)
                        {
                            for (; row < row_end; row++) {
                                for (size_t col = 0; col < ncols_; col++) {
                                    if (std::isnan((float)X[(size_t)row + col*nrows_]) ||
                                        std::isinf((float)X[(size_t)row + col*nrows_])
                                    ) {
                                        indices_.push_back(col);
                                        values_.push_back(
                                            NAN
                                        );
                                    }
                                }
                                indptr[row+1] = indices_.size();
                            }
                        }

                        else
                        {
                            for (; row < row_end; row++) {
                                for (size_t col = 0; col < ncols_; col++) {
                                    if (X[(size_t)row + col*nrows_] == NA_INTEGER
                                    ) {
                                        indices_.push_back(col);
                                        values_.push_back(
                                            std::isinf((double)X[(size_t)row + col*nrows_])?
                                                NAN : X[(size_t)row + col*nrows_]
                                        );
                                    }
                                }
                                indptr[row+1] = indices_.size();
                            }
                        }

                        break;

                    }

                    else if (row == *curr_i-1 + offset) {

                        if (simple_mult)
                        {
                            for (size_t col = 0; col < ncols_; col++) {
                                indices_.push_back(col);
                                values_.push_back(X[(size_t)row + col*nrows_] * xx[curr_i - ii_begin]);
                            }
                        }

                        else
                        {
                            for (size_t col = 0; col < ncols_; col++) {
                                indices_.push_back(col);
                                values_.push_back(
                                    (X[(size_t)row + col*nrows_] == NA_INTEGER)? NAN : 
                                    (X[(size_t)row + col*nrows_] * xx[curr_i - ii_begin])
                                );
                            }
                        }
                        indptr[row+1] = indices_.size();
                        row++;
                        curr_i++;

                    }

                    else if (row > *curr_i-1 + offset) {

                        curr_i = std::lower_bound(curr_i, end_i, row+1 - offset);

                    }

                    else {

                        if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value)
                        {
                            for (; row < *curr_i-1 + offset; row++) {
                                for (size_t col = 0; col < ncols_; col++) {
                                    if (ISNAN(X[(size_t)row + col*nrows_]) ||
                                        std::isinf((double)X[(size_t)row + col*nrows_])
                                    ) {
                                        indices_.push_back(col);
                                        values_.push_back(
                                            std::isinf((double)X[(size_t)row + col*nrows_])?
                                                NAN : X[(size_t)row + col*nrows_]
                                        );
                                    }
                                }
                                indptr[row+1] = indices_.size();
                            }
                        }

                        else if (std::is_same<InputDType, float>::value)
                        {
                            for (; row < *curr_i-1 + offset; row++) {
                                for (size_t col = 0; col < ncols_; col++) {
                                    if (std::isnan((float)X[(size_t)row + col*nrows_]) ||
                                        std::isinf((float)X[(size_t)row + col*nrows_])
                                    ) {
                                        indices_.push_back(col);
                                        values_.push_back(
                                            NAN
                                        );
                                    }
                                }
                                indptr[row+1] = indices_.size();
                            }
                        }

                        else
                        {
                            for (; row < *curr_i-1 + offset; row++) {
                                for (size_t col = 0; col < ncols_; col++) {
                                    if (X[(size_t)row + col*nrows_] == NA_INTEGER
                                    ) {
                                        indices_.push_back(col);
                                        values_.push_back(
                                            NA_REAL
                                        );
                                    }
                                }
                                indptr[row+1] = indices_.size();
                            }
                        }

                    }
                }
            }


            VectorConstructorArgs args;
            args.from_cpp_vec = true; args.as_integer = true; args.int_vec_from = &indices_;
            indices = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
            indices_.clear(); indices_.shrink_to_fit();
            args.as_integer = false; args.num_vec_from = &values_;
            values = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
        }

        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }

    else
    {
        size_t ix_max = (size_t)nrows * (size_t)ncols;
        size_t n_recycles = ix_max / length + ((ix_max % length) != 0);

        size_t ix_curr;
        Rcpp::NumericMatrix out_(nrows, ncols);
        auto out = out_.begin();
        double val;

        if (keep_NAs)
        {
            size_t tot = (size_t)nrows * (size_t)ncols;

            if (std::is_same<RcppMatrix, Rcpp::NumericMatrix>::value) {
                for (size_t ix = 0; ix < tot; ix++)
                    out[ix] = ISNAN(X[ix])? X[ix] : (std::isinf((double)X[ix])? NAN : 0);
            }

            else if (std::is_same<InputDType, float>::value) {
                for (size_t ix = 0; ix < tot; ix++)
                    out[ix] = (std::isnan((float)X[ix]) || std::isinf((float)X[ix])? NAN : 0);
            }

            else {
                for (size_t ix = 0; ix < tot; ix++)
                    out[ix] = (X[ix] == NA_INTEGER)? NA_REAL : 0;
            }
        }

        if (simple_mult)
        {
            for (size_t ix = 0; ix < nnz; ix++)
            {
                val = xx[ix];
                for (size_t recycle = 0; recycle < n_recycles; recycle++)
                {
                    ix_curr = ii[ix] - 1 + recycle*length;
                    if (ix_curr > ix_max)
                        break;
                    row = ix_curr % nrows_;
                    col = ix_curr / nrows_;
                    out[row + col*nrows_] = X[row + col*nrows_] * val;
                }
            }
        }

        else
        {
            for (size_t ix = 0; ix < nnz; ix++)
            {
                val = xx[ix];
                for (size_t recycle = 0; recycle < n_recycles; recycle++)
                {
                    ix_curr = ii[ix] - 1 + recycle*length;
                    if (ix_curr > ix_max)
                        break;
                    row = ix_curr % nrows_;
                    col = ix_curr / nrows_;
                    out[row + col*nrows_] = (X[row + col*nrows_] == NA_INTEGER)? NAN : (X[row + col*nrows_] * val);
                }
            }
        }

        return Rcpp::List::create(Rcpp::_["X_dense"] = out_);
    }

    return Rcpp::List(); /* <- won't be reached, but CRAN may complain otherwise */
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_elemwise_dense_by_svec_numeric
(
    Rcpp::NumericMatrix X_,
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const int length,
    const int keep_NAs
)
{
    return multiply_elemwise_dense_by_svec_template<Rcpp::NumericMatrix, double>(
        X_,
        ii,
        xx,
        length,
        keep_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_elemwise_dense_by_svec_integer
(
    Rcpp::IntegerMatrix X_,
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const int length,
    const int keep_NAs
)
{
    return multiply_elemwise_dense_by_svec_template<Rcpp::IntegerMatrix, int>(
        X_,
        ii,
        xx,
        length,
        keep_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_elemwise_dense_by_svec_logical
(
    Rcpp::LogicalMatrix X_,
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const int length,
    const int keep_NAs
)
{
    return multiply_elemwise_dense_by_svec_template<Rcpp::LogicalMatrix, int>(
        X_,
        ii,
        xx,
        length,
        keep_NAs
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List multiply_elemwise_dense_by_svec_float32
(
    Rcpp::IntegerMatrix X_,
    Rcpp::IntegerVector ii,
    Rcpp::NumericVector xx,
    const int length,
    const int keep_NAs
)
{
    return multiply_elemwise_dense_by_svec_template<Rcpp::IntegerMatrix, float>(
        X_,
        ii,
        xx,
        length,
        keep_NAs
    );
}

#ifdef __clang__
#   pragma clang diagnostic pop
#endif
