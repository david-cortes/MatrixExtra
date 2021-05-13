#include "MatrixExtra.h"

#ifdef __clang__
#   pragma clang diagnostic push
#   pragma clang diagnostic ignored "-Wpass-failed=transform-warning"
#endif

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector concat_indptr2(Rcpp::IntegerVector ptr1, Rcpp::IntegerVector ptr2)
{
    Rcpp::IntegerVector out(ptr1.size() + ptr2.size() - 1);
    std::copy(INTEGER(ptr1), INTEGER(ptr1) + ptr1.size(), INTEGER(out));
    size_t st_second = ptr1.size();
    int offset = ptr1[ptr1.size()-1];
    #ifdef _OPENMP
    #pragma omp simd
    #endif
    for (size_t row = 1; row < (size_t)ptr2.size(); row++)
        out[st_second + row - 1] = offset + ptr2[row];
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::S4 concat_csr_batch(Rcpp::ListOf<Rcpp::S4> objects, Rcpp::S4 out)
{
    size_t n_inputs = objects.size();
    RbindedType otype;
    if (out.inherits("ngRMatrix"))
        otype = ngRMatrix;
    else if (out.inherits("lgRMatrix"))
        otype = lgRMatrix;
    else
        otype = dgRMatrix;

    int *indptr_out = INTEGER(out.slot("p"));
    int *indices_out = INTEGER(out.slot("j"));
    double *values_out = (otype == dgRMatrix)? REAL(out.slot("x")) : nullptr;
    int *values_out_bool = (otype == lgRMatrix)? LOGICAL(out.slot("x")) : nullptr;

    int nrows_add, nnz_add;
    int *indptr_obj, *indices_obj;
    double *xvals_obj = nullptr;
    int *lvals_obj = nullptr;
    int *ivals_obj = nullptr;

    indptr_out[0] = 0;
    int curr_pos = 0;
    int curr_row = 0;

    for (size_t ix = 0; ix < n_inputs; ix++)
    {

        if (objects[ix].hasSlot("j"))
        {
            
            indptr_obj = INTEGER(objects[ix].slot("p"));
            indices_obj = INTEGER(objects[ix].slot("j"));
            nnz_add = Rf_xlength(objects[ix].slot("j"));
            xvals_obj = objects[ix].inherits("dgRMatrix")? REAL(objects[ix].slot("x")) : nullptr;
            lvals_obj = objects[ix].inherits("lgRMatrix")? LOGICAL(objects[ix].slot("x")) : nullptr;
            nrows_add = INTEGER(objects[ix].slot("Dim"))[0];

            /* indptr */
            for (int row = 0; row < nrows_add; row++)
                indptr_out[row + curr_row + 1] = indptr_out[curr_row] + indptr_obj[row+1];
            curr_row += nrows_add;
            /* indices */
            std::copy(indices_obj, indices_obj + nnz_add, indices_out + curr_pos);
            /* values, if applicable */
            if (otype == dgRMatrix) {
                if (xvals_obj != nullptr)
                    std::copy(xvals_obj, xvals_obj + nnz_add, values_out + curr_pos);
                else if (lvals_obj != nullptr)
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int el = 0; el < nnz_add; el++)
                        values_out[el + curr_pos] = (lvals_obj[el] == NA_LOGICAL)? (NA_REAL) : lvals_obj[el];
                else
                    std::fill(values_out + curr_pos, values_out + curr_pos + nnz_add, 1.);
            }
            
            else if (otype == lgRMatrix) {
                if (lvals_obj != nullptr)
                    std::copy(lvals_obj, lvals_obj + nnz_add, values_out_bool + curr_pos);
                else if (xvals_obj != nullptr)
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int el = 0; el < nnz_add; el++)
                        values_out_bool[el + curr_pos] = ISNAN(xvals_obj[el])? NA_LOGICAL : (bool)xvals_obj[el];
                else
                    std::fill(values_out_bool + curr_pos, values_out_bool + curr_pos + nnz_add, (int)true);
            }
        }

        else
        {
            
            indices_obj = INTEGER(objects[ix].slot("i"));
            nnz_add = Rf_xlength(objects[ix].slot("i"));
            indptr_out[curr_row + 1] = indptr_out[curr_row] + nnz_add;
            curr_row++;
            #ifdef _OPENMP
            #pragma omp simd
            #endif
            for (int el = 0; el < nnz_add; el++)
                indices_out[el + curr_pos] = indices_obj[el] - 1;

            if (otype == dgRMatrix) {
                
                if (objects[ix].inherits("dsparseVector")) {
                    xvals_obj = REAL(objects[ix].slot("x"));
                    std::copy(xvals_obj, xvals_obj + nnz_add, values_out + curr_pos);
                } else if (objects[ix].inherits("isparseVector")) {
                    ivals_obj = INTEGER(objects[ix].slot("x"));
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int el = 0; el < nnz_add; el++)
                        values_out[el + curr_pos] = (ivals_obj[el] == NA_INTEGER)? NA_REAL : ivals_obj[el];
                } else if (objects[ix].inherits("lsparseVector")) {
                    lvals_obj = LOGICAL(objects[ix].slot("x"));
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int el = 0; el < nnz_add; el++)
                        values_out[el + curr_pos] = (lvals_obj[el] == NA_LOGICAL)? NA_REAL : (bool)lvals_obj[el];
                } else if (objects[ix].inherits("nsparseVector")) {
                    std::fill(values_out + curr_pos, values_out + curr_pos + nnz_add, 1.);
                } else {
                    char errmsg[100];
                    std::snprintf(errmsg, 99, "Invalid vector type in argument %d.\n", (int)ix);
                    Rcpp::stop(errmsg);
                }

            }

            else if (otype == lgRMatrix) {

                if (objects[ix].inherits("dsparseVector")) {
                    xvals_obj = REAL(objects[ix].slot("x"));
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int el = 0; el < nnz_add; el++)
                        values_out_bool[el + curr_pos] = ISNAN(xvals_obj[el])? NA_LOGICAL : (bool)xvals_obj[el];
                } else if (objects[ix].inherits("isparseVector")) {
                    ivals_obj = INTEGER(objects[ix].slot("x"));
                    #ifdef _OPENMP
                    #pragma omp simd
                    #endif
                    for (int el = 0; el < nnz_add; el++)
                        values_out_bool[el + curr_pos] = (ivals_obj[el] == NA_INTEGER)? NA_LOGICAL : (bool)ivals_obj[el];
                } else if (objects[ix].inherits("lsparseVector")) {
                    lvals_obj = LOGICAL(objects[ix].slot("x"));
                    std::copy(lvals_obj, lvals_obj + nnz_add, values_out_bool + curr_pos);
                } else if (objects[ix].inherits("nsparseVector")) {
                    std::fill(values_out_bool + curr_pos, values_out_bool + curr_pos + nnz_add, (int)true);
                } else {
                    char errmsg[100];
                    std::snprintf(errmsg, 99, "Invalid vector type in argument %d.\n", (int)ix);
                    Rcpp::stop(errmsg);
                }

            }
        }

        curr_pos += nnz_add;
    }

    return out;
}

#ifdef __clang__
#   pragma clang diagnostic pop
#endif
