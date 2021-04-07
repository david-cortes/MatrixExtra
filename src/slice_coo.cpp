#include "MatrixExtra.h"

// [[Rcpp::export(rng = false)]]
int slice_coo_single
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    int i, int j
)
{
    size_t nnz = ii.size();
    for (size_t ix = 0; ix < nnz; ix++)
    {
        if (ii[ix] == i && jj[ix] == j)
        {
            return ix + 1;
        }
    }
    return 0;
}

void process_i_arbitrary
(
    Rcpp::IntegerVector rows_take_base1,
    const bool all_i, const bool i_is_seq, const bool i_is_rev_seq,
    const int first_i, const int last_i,
    const int nrows,
    int &min_i, int &max_i,
    std::unordered_map<int, int> &i_mapping,
    std::unordered_map<int, std::vector<int>> &i_indices_rep,
    bool &i_has_duplicates
)
{
    if (all_i) {
        min_i = 0;
        max_i = nrows - 1;
    } else if (i_is_seq) {
        min_i = first_i;
        max_i = last_i;
    } else if (i_is_rev_seq) {
        min_i = last_i;
        max_i = first_i;
    } else {
        min_i = *std::min_element(rows_take_base1.begin(), rows_take_base1.end()) - 1;
        max_i = *std::max_element(rows_take_base1.begin(), rows_take_base1.end()) - 1;
    }


    if (!all_i && !i_is_seq && !i_is_rev_seq)
    {
        for (int ix = 0; ix < rows_take_base1.size(); ix++)
            i_mapping[rows_take_base1[ix]-1] = ix;
        i_has_duplicates = i_mapping.size() != rows_take_base1.size();

        if (i_has_duplicates)
        {
            for (int ix = 0; ix < rows_take_base1.size(); ix++)
                i_indices_rep[rows_take_base1[ix]-1].push_back(ix);
        }        
    }

    else {
        i_has_duplicates = false;
    }
}

void post_process_seq
(
    int *restrict ii_out, size_t n,
    const bool all_i, const bool i_is_seq, const bool i_is_rev_seq,
    const int first_i, const int last_i
)
{
    if (!all_i)
    {
        if (i_is_seq) {
            #pragma omp simd
            for (size_t ix = 0; ix < n; ix++)
                ii_out[ix] -= first_i;
        }
        else if (i_is_rev_seq) {
            #pragma omp simd
            for (size_t ix = 0; ix < n; ix++)
                ii_out[ix] = first_i - ii_out[ix];
        }
    }
}

template <class RcppVector, class InputDType, class CompileFlag>
Rcpp::List slice_coo_arbitrary_template
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    RcppVector xx,
    Rcpp::IntegerVector rows_take_base1,
    Rcpp::IntegerVector cols_take_base1,
    bool all_i, bool all_j,
    bool i_is_seq, bool j_is_seq,
    bool i_is_rev_seq, bool j_is_rev_seq,
    int nrows, int ncols
)
{
    size_t nnz = ii.size();
    size_t size_reserve = get_size_reserve(nnz, rows_take_base1.size(), cols_take_base1.size());


    size_t curr = 0;
    int first_i = rows_take_base1[0] - 1;
    int first_j = cols_take_base1[0] - 1;
    int last_i  = rows_take_base1[rows_take_base1.size()-1] - 1;
    int last_j  = cols_take_base1[cols_take_base1.size()-1] - 1;
    std::vector<size_t> take;


    if ((all_i || i_is_seq) && (all_j || j_is_seq))
    {
        take = std::vector<size_t>(size_reserve);
        int min_i = std::min(first_i, last_i), max_i = std::max(first_i, last_i);
        int min_j = std::min(first_j, last_j), max_j = std::max(first_j, last_j);

        if (!all_i && !all_j) {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
                    take[curr++] = ix;
            }
        }

        else if (all_i) {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (jj[ix] >= min_j && jj[ix] <= max_j)
                    take[curr++] = ix;
            }
        }

        else {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (ii[ix] >= min_i && ii[ix] <= max_i)
                    take[curr++] = ix;
            }
        }

        goto easy_end;
    }

    else if ((all_i || i_is_seq || i_is_rev_seq) && (all_j || j_is_seq || j_is_rev_seq))
    {
        take = std::vector<size_t>(size_reserve);
        int min_i = std::min(first_i, last_i), max_i = std::max(first_i, last_i);
        int min_j = std::min(first_j, last_j), max_j = std::max(first_j, last_j);

        if (all_i) {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (jj[ix] >= min_j && jj[ix] <= max_j)
                    take[curr++] = ix;
            }
        }

        else if (all_j) {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (ii[ix] >= min_i && ii[ix] <= max_i)
                    take[curr++] = ix;
            }
        }

        else if (i_is_seq) {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
                    take[curr++] = ix;
            }
        }

        else if (j_is_seq) {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
                    take[curr++] = ix;
            }
        }

        else {
            for (size_t ix = 0; ix < nnz; ix++) {
                if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
                    take[curr++] = ix;
            }
        }

        goto easy_end;
    }

    if (false)
    {
        easy_end:
        Rcpp::IntegerVector ii_out(curr);
        Rcpp::IntegerVector jj_out(curr);
        RcppVector xx_out(std::is_same<CompileFlag, bool>::value? curr : 0);

        for (size_t ix = 0; ix < curr; ix++) ii_out[ix] = ii[take[ix]];
        for (size_t ix = 0; ix < curr; ix++) jj_out[ix] = jj[take[ix]];
        if (std::is_same<CompileFlag, bool>::value)
        for (size_t ix = 0; ix < curr; ix++) xx_out[ix] = xx[take[ix]];

        post_process_seq(
            INTEGER(ii_out), curr,
            all_i, i_is_seq, i_is_rev_seq,
            first_i, last_i
        );
        post_process_seq(
            INTEGER(jj_out), curr,
            all_j, j_is_seq, j_is_rev_seq,
            first_j, last_j
        );

        return Rcpp::List::create(
            Rcpp::_["ii"] = ii_out,
            Rcpp::_["jj"] = jj_out,
            Rcpp::_["xx"] = xx_out
        );
    }




    int min_i, max_i, min_j, max_j;
    std::unordered_map<int, int> i_mapping, j_mapping;
    std::unordered_map<int, std::vector<int>> i_indices_rep, j_indices_rep;
    bool i_has_duplicates = false, j_has_duplicates = false;

    process_i_arbitrary(
        rows_take_base1,
        all_i, i_is_seq, i_is_rev_seq,
        first_i, last_i,
        nrows,
        min_i, max_i,
        i_mapping,
        i_indices_rep,
        i_has_duplicates
    );

    process_i_arbitrary(
        cols_take_base1,
        all_j, j_is_seq, j_is_rev_seq,
        first_j, last_j,
        ncols,
        min_j, max_j,
        j_mapping,
        j_indices_rep,
        j_has_duplicates
    );

    curr = 0;
    std::vector<int> ii_out;
    std::vector<int> jj_out;
    std::vector<InputDType> xx_out;

    if (!i_has_duplicates && !j_has_duplicates) {
        ii_out = std::vector<int>(size_reserve);
        jj_out = std::vector<int>(size_reserve);
        if (std::is_same<CompileFlag, bool>::value)
        xx_out = std::vector<InputDType>(size_reserve);
    }

    else {
        ii_out.reserve(size_reserve);
        jj_out.reserve(size_reserve);
        if (std::is_same<CompileFlag, bool>::value)
        xx_out.reserve(size_reserve);
    }

    std::unordered_map<int, int>::iterator res_i, res_j;
    std::unordered_map<int, std::vector<int>>::iterator res_in, res_jn;
    bool post_process = false;

    if ((all_i || i_is_seq || i_is_rev_seq) && !j_has_duplicates) {

        for (size_t ix = 0; ix < nnz; ix++)
        {
            if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
            {
                res_j = j_mapping.find(jj[ix]);
                if (res_j != j_mapping.end())
                {
                    ii_out[curr] = ii[ix];
                    jj_out[curr] = res_j->second;
                    if (std::is_same<CompileFlag, bool>::value)
                    xx_out[curr] = xx[ix];
                    curr++;
                }
            }
        }

        post_process = true;

    }

    else if ((all_j || j_is_seq || j_is_rev_seq) && !i_has_duplicates) {

        for (size_t ix = 0; ix < nnz; ix++)
        {
            if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
            {
                res_i = i_mapping.find(ii[ix]);
                if (res_i != i_mapping.end())
                {
                    ii_out[curr] = res_i->second;
                    jj_out[curr] = jj[ix];
                    if (std::is_same<CompileFlag, bool>::value)
                    xx_out[curr] = xx[ix];
                    curr++;
                }
            }
        }

        post_process = true;

    }

    else if (!i_has_duplicates && !j_has_duplicates) {
        
        for (size_t ix = 0; ix < nnz; ix++)
        {
            if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
            {
                res_i = i_mapping.find(ii[ix]);
                if (res_i != i_mapping.end())
                {
                    res_j = j_mapping.find(jj[ix]);
                    if (res_j != j_mapping.end())
                    {
                        ii_out[curr] = res_i->second;
                        jj_out[curr] = res_j->second;
                        if (std::is_same<CompileFlag, bool>::value)
                        xx_out[curr] = xx[ix];
                        curr++;
                    }
                }
            }
        }

    }

    else if (!i_has_duplicates) {

        for (size_t ix = 0; ix < nnz; ix++)
        {
            if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
            {
                res_i = i_mapping.find(ii[ix]);
                if (res_i != i_mapping.end())
                {
                    res_jn = j_indices_rep.find(jj[ix]);
                    if (res_jn != j_indices_rep.end())
                    {
                        for (auto el : res_jn->second)
                        {
                            ii_out.push_back(res_i->second);
                            jj_out.push_back(el);
                            if (std::is_same<CompileFlag, bool>::value)
                            xx_out.push_back(xx[ix]);
                        }
                    }
                }
            }
        }

    }

    else if (!j_has_duplicates) {

        for (size_t ix = 0; ix < nnz; ix++)
        {
            if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
            {
                res_j = j_mapping.find(jj[ix]);
                if (res_j != j_mapping.end())
                {
                    res_in = i_indices_rep.find(ii[ix]);
                    if (res_in != i_indices_rep.end())
                    {
                        for (auto el : res_in->second)
                        {
                            ii_out.push_back(el);
                            jj_out.push_back(res_j->second);
                            if (std::is_same<CompileFlag, bool>::value)
                            xx_out.push_back(xx[ix]);
                        }
                    }
                }
            }
        }

    }

    else {

        for (size_t ix = 0; ix < nnz; ix++)
        {
            if (ii[ix] >= min_i && ii[ix] <= max_i && jj[ix] >= min_j && jj[ix] <= max_j)
            {
                res_in = i_indices_rep.find(ii[ix]);
                if (res_in != i_indices_rep.end())
                {
                    res_jn = j_indices_rep.find(jj[ix]);
                    if (res_jn != j_indices_rep.end())
                    {
                        for (auto el_i : res_in->second)
                        {
                            for (auto el_j : res_jn->second)
                            {
                                ii_out.push_back(el_i);
                                jj_out.push_back(el_j);
                                if (std::is_same<CompileFlag, bool>::value)
                                xx_out.push_back(xx[ix]);
                            }
                        }
                    }
                }
            }
        }

    }


    if (post_process)
    {
        post_process_seq(
            ii_out.data(), curr,
            all_i, i_is_seq, i_is_rev_seq,
            first_i, last_i
        );
        post_process_seq(
            jj_out.data(), curr,
            all_j, j_is_seq, j_is_rev_seq,
            first_j, last_j
        );
    }

    if (i_has_duplicates || j_has_duplicates)
        curr = ii_out.size();

    Rcpp::IntegerVector ii_out_(ii_out.begin(), ii_out.begin() + curr);
    Rcpp::IntegerVector jj_out_(jj_out.begin(), jj_out.begin() + curr);
    RcppVector xx_out_;
    if (std::is_same<CompileFlag, bool>::value)
    xx_out_ = RcppVector(xx_out.begin(), xx_out.begin() + curr);


    return Rcpp::List::create(
        Rcpp::_["ii"] = ii_out_,
        Rcpp::_["jj"] = jj_out_,
        Rcpp::_["xx"] = xx_out_
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List slice_coo_arbitrary_numeric
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::NumericVector xx,
    Rcpp::IntegerVector rows_take_base1,
    Rcpp::IntegerVector cols_take_base1,
    bool all_i, bool all_j,
    bool i_is_seq, bool j_is_seq,
    bool i_is_rev_seq, bool j_is_rev_seq,
    int nrows, int ncols
)
{
    return slice_coo_arbitrary_template<Rcpp::NumericVector, double, bool>(
        ii,
        jj,
        xx,
        rows_take_base1,
        cols_take_base1,
        all_i, all_j,
        i_is_seq, j_is_seq,
        i_is_rev_seq, j_is_rev_seq,
        nrows, ncols
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List slice_coo_arbitrary_logical
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::LogicalVector xx,
    Rcpp::IntegerVector rows_take_base1,
    Rcpp::IntegerVector cols_take_base1,
    bool all_i, bool all_j,
    bool i_is_seq, bool j_is_seq,
    bool i_is_rev_seq, bool j_is_rev_seq,
    int nrows, int ncols
)
{
    return slice_coo_arbitrary_template<Rcpp::LogicalVector, int, bool>(
        ii,
        jj,
        xx,
        rows_take_base1,
        cols_take_base1,
        all_i, all_j,
        i_is_seq, j_is_seq,
        i_is_rev_seq, j_is_rev_seq,
        nrows, ncols
    );
}

// [[Rcpp::export(rng = false)]]
Rcpp::List slice_coo_arbitrary_binary
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::IntegerVector rows_take_base1,
    Rcpp::IntegerVector cols_take_base1,
    bool all_i, bool all_j,
    bool i_is_seq, bool j_is_seq,
    bool i_is_rev_seq, bool j_is_rev_seq,
    int nrows, int ncols
)
{
    return slice_coo_arbitrary_template<Rcpp::NumericVector, double, int>(
        ii,
        jj,
        Rcpp::NumericVector(),
        rows_take_base1,
        cols_take_base1,
        all_i, all_j,
        i_is_seq, j_is_seq,
        i_is_rev_seq, j_is_rev_seq,
        nrows, ncols
    );
}
