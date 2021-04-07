#include "MatrixExtra.h"

template <class RcppVector>
Rcpp::List cbind_csr
(
    const Rcpp::IntegerVector X_csr_indptr,
    const Rcpp::IntegerVector X_csr_indices,
    const RcppVector X_csr_values,
    const Rcpp::IntegerVector Y_csr_indptr,
    const Rcpp::IntegerVector Y_csr_indices_plus_ncol,
    const RcppVector Y_csr_values
)
{
    int nrows = std::max(X_csr_indptr.size(), Y_csr_indptr.size()) - 1;
    int nrows_min = std::min(X_csr_indptr.size(), Y_csr_indptr.size()) - 1;
    Rcpp::IntegerVector indptr(nrows + 1);
    Rcpp::IntegerVector indices(X_csr_indices.size() + Y_csr_indices_plus_ncol.size());
    RcppVector values;
    if (X_csr_values.size() || Y_csr_values.size())
        values = RcppVector(indices.size());

    if (!indices.size())
    {
        return Rcpp::List::create(
            Rcpp::_["indptr"] = indptr,
            Rcpp::_["indices"] = indices,
            Rcpp::_["values"] = values
        );
    }
    
    for (int row = 0; row < nrows_min; row++)
        indptr[row+1] = indptr[row] + X_csr_indptr[row+1] - X_csr_indptr[row] + Y_csr_indptr[row+1] - Y_csr_indptr[row];
    if (X_csr_indptr.size() > Y_csr_indptr.size()) {
        for (int row = nrows_min; row < nrows; row++)
            indptr[row+1] = indptr[row] + X_csr_indptr[row+1] - X_csr_indptr[row];
    } else if (Y_csr_indptr.size() > X_csr_indptr.size()) {
        for (int row = nrows_min; row < nrows; row++)
            indptr[row+1] = indptr[row] + Y_csr_indptr[row+1] - Y_csr_indptr[row];
    }

    const bool has_values = values.size() > 0;
    int nnz_x;
    const auto X_csr_indices_begin = X_csr_indices.begin();
    const auto Y_csr_indices_plus_ncol_begin = Y_csr_indices_plus_ncol.begin();
    const auto X_csr_values_begin = X_csr_values.begin();
    const auto Y_csr_values_begin = Y_csr_values.begin();
    const auto indices_begin = indices.begin();
    const auto values_begin = values.begin();

    /* TODO: check if this runs faster when parallelized */

    for (int row = 0; row < nrows_min; row++)
    {
        nnz_x = X_csr_indptr[row+1] - X_csr_indptr[row];
        std::copy(X_csr_indices_begin + X_csr_indptr[row],
                  X_csr_indices_begin + X_csr_indptr[row+1],
                  indices_begin + indptr[row]);
        std::copy(Y_csr_indices_plus_ncol_begin + Y_csr_indptr[row],
                  Y_csr_indices_plus_ncol_begin + Y_csr_indptr[row+1],
                  indices_begin + indptr[row] + nnz_x);

        if (has_values)
        {
        std::copy(X_csr_values_begin + X_csr_indptr[row],
                  X_csr_values_begin + X_csr_indptr[row+1],
                  values_begin + indptr[row]);
        std::copy(Y_csr_values_begin + Y_csr_indptr[row],
                  Y_csr_values_begin + Y_csr_indptr[row+1],
                  values_begin + indptr[row] + nnz_x);
        }
    }

    if (X_csr_indptr.size() > Y_csr_indptr.size())
    {
        std::copy(X_csr_indices_begin + X_csr_indptr[nrows_min],
                  X_csr_indices.end(),
                  indices_begin + indptr[nrows_min]);
        if (has_values)
        std::copy(X_csr_values_begin + X_csr_indptr[nrows_min],
                  X_csr_values.end(),
                  values_begin + indptr[nrows_min]);
    }

    else if (Y_csr_indptr.size() > X_csr_indptr.size()) {
        std::copy(Y_csr_indices_plus_ncol_begin + Y_csr_indptr[nrows_min],
                  Y_csr_indices_plus_ncol.end(),
                  indices_begin + indptr[nrows_min]);
        if (has_values)
        std::copy(Y_csr_values_begin + Y_csr_indptr[nrows_min],
                  Y_csr_values.end(),
                  values_begin + indptr[nrows_min]);
    }

    return Rcpp::List::create(
        Rcpp::_["indptr"] = indptr,
        Rcpp::_["indices"] = indices,
        Rcpp::_["values"] = values
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cbind_csr_numeric
(
    Rcpp::IntegerVector X_csr_indptr,
    Rcpp::IntegerVector X_csr_indices,
    Rcpp::NumericVector X_csr_values,
    Rcpp::IntegerVector Y_csr_indptr,
    Rcpp::IntegerVector Y_csr_indices_plus_ncol,
    Rcpp::NumericVector Y_csr_values
)
{
    return cbind_csr(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        Y_csr_indptr,
        Y_csr_indices_plus_ncol,
        Y_csr_values
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cbind_csr_logical
(
    Rcpp::IntegerVector X_csr_indptr,
    Rcpp::IntegerVector X_csr_indices,
    Rcpp::LogicalVector X_csr_values,
    Rcpp::IntegerVector Y_csr_indptr,
    Rcpp::IntegerVector Y_csr_indices_plus_ncol,
    Rcpp::LogicalVector Y_csr_values
)
{
    return cbind_csr(
        X_csr_indptr,
        X_csr_indices,
        X_csr_values,
        Y_csr_indptr,
        Y_csr_indices_plus_ncol,
        Y_csr_values
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cbind_csr_binary
(
    Rcpp::IntegerVector X_csr_indptr,
    Rcpp::IntegerVector X_csr_indices,
    Rcpp::IntegerVector Y_csr_indptr,
    Rcpp::IntegerVector Y_csr_indices_plus_ncol
)
{
    return cbind_csr(
        X_csr_indptr,
        X_csr_indices,
        Rcpp::NumericVector(),
        Y_csr_indptr,
        Y_csr_indices_plus_ncol,
        Rcpp::NumericVector()
    );
}
