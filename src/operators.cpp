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
