#include "MatrixExtra.h"

template <class RcppVector, class InputDType, class CompileFlag>
InputDType slice_coo_single_template
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    RcppVector xx,
    const int i, const int j
)
{
    const size_t nnz = ii.size();
    for (size_t ix = 0; ix < nnz; ix++)
    {
        if (ii[ix] == i && jj[ix] == j)
        {
            return std::is_same<CompileFlag, bool>::value? xx[ix] : 1;
        }
    }
    return 0;
}

// [[Rcpp::export(rng = false)]]
double slice_coo_single_numeric
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::NumericVector xx,
    int i, int j
)
{
    return slice_coo_single_template<Rcpp::NumericVector, double, bool>(
        ii,
        jj,
        xx,
        i, j
    );
}

// [[Rcpp::export(rng = false)]]
bool slice_coo_single_logical
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    Rcpp::LogicalVector xx,
    int i, int j
)
{
    return slice_coo_single_template<Rcpp::LogicalVector, int, bool>(
        ii,
        jj,
        xx,
        i, j
    );
}

// [[Rcpp::export(rng = false)]]
bool slice_coo_single_binary
(
    Rcpp::IntegerVector ii,
    Rcpp::IntegerVector jj,
    int i, int j
)
{
    return slice_coo_single_template<bool*, bool, int>(
        ii,
        jj,
        nullptr,
        i, j
    );
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
    std::unique_ptr<size_t[]> take;


    if ((all_i || i_is_seq) && (all_j || j_is_seq))
    {
        take = std::unique_ptr<size_t[]>(new size_t[size_reserve]);
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
        take = std::unique_ptr<size_t[]>(new size_t[size_reserve]);
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
        VectorConstructorArgs args;
        args.as_integer = true; args.from_cpp_vec = false; args.as_logical = false; args.size = curr;
        Rcpp::IntegerVector ii_out = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
        Rcpp::IntegerVector jj_out = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
        RcppVector xx_out;
        if (std::is_same<CompileFlag, bool>::value)
        {
            if (std::is_same<RcppVector, Rcpp::LogicalVector>::value) {
                args.as_logical = true; args.as_integer = true;
            }
            else if (std::is_same<RcppVector, Rcpp::IntegerVector>::value) {
                args.as_logical = false; args.as_integer = true;
            } else {
                args.as_integer = false;
            }

            xx_out = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
        }

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
    std::unique_ptr<int[]> ii_out_;
    std::unique_ptr<int[]> jj_out_;
    std::unique_ptr<InputDType[]> xx_out_;
    bool use_vectors = false;

    if (!i_has_duplicates && !j_has_duplicates) {
        ii_out = std::vector<int>(size_reserve);
        jj_out = std::vector<int>(size_reserve);
        if (std::is_same<CompileFlag, bool>::value)
        xx_out = std::vector<InputDType>(size_reserve);

        use_vectors = false;
        ii_out_ = std::unique_ptr<int[]>(new int[size_reserve]);
        jj_out_ = std::unique_ptr<int[]>(new int[size_reserve]);
        if (std::is_same<CompileFlag, bool>::value)
        xx_out_ = std::unique_ptr<InputDType[]>(new InputDType[size_reserve]);
    }

    else {
        use_vectors = true;
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
                    ii_out_[curr] = ii[ix];
                    jj_out_[curr] = res_j->second;
                    if (std::is_same<CompileFlag, bool>::value)
                    xx_out_[curr] = xx[ix];
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
                    ii_out_[curr] = res_i->second;
                    jj_out_[curr] = jj[ix];
                    if (std::is_same<CompileFlag, bool>::value)
                    xx_out_[curr] = xx[ix];
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
                        ii_out_[curr] = res_i->second;
                        jj_out_[curr] = res_j->second;
                        if (std::is_same<CompileFlag, bool>::value)
                        xx_out_[curr] = xx[ix];
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
            use_vectors? ii_out.data() : ii_out_.get(), curr,
            all_i, i_is_seq, i_is_rev_seq,
            first_i, last_i
        );
        post_process_seq(
            use_vectors? jj_out.data() : jj_out_.get(), curr,
            all_j, j_is_seq, j_is_rev_seq,
            first_j, last_j
        );
    }

    if (i_has_duplicates || j_has_duplicates)
        curr = ii_out.size();

    VectorConstructorArgs args;
    args.as_integer = true; args.as_logical = false;
    args.size = curr; args.cpp_lim_size = true;

    if (use_vectors) {
        args.from_cpp_vec = true;
        args.int_vec_from = &ii_out;
    } else {
        args.from_cpp_vec = false;
        args.from_pointer = true;
        args.int_pointer_from = ii_out_.get();
    }
    Rcpp::IntegerVector ii_out__ = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    ii_out.clear(); ii_out.shrink_to_fit(); ii_out_.reset();
    args.int_vec_from = &jj_out; args.int_pointer_from = jj_out_.get();
    Rcpp::IntegerVector jj_out__ = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    jj_out.clear(); jj_out.shrink_to_fit(); jj_out_.reset();
    RcppVector xx_out__;
    if (std::is_same<CompileFlag, bool>::value)
    {
        if (std::is_same<RcppVector, Rcpp::LogicalVector>::value) {
            args.as_integer = true; args.as_logical = true;
            if (use_vectors)
                args.int_vec_from = &xx_out;
            else
                args.int_pointer_from = xx_out_.get();
        }
        else if (std::is_same<RcppVector, Rcpp::IntegerVector>::value) {
            args.as_integer = true; args.as_logical = false;
            if (use_vectors)
                args.int_vec_from = &xx_out;
            else
                args.int_pointer_from = xx_out_.get();
        }
        else {
            args.as_integer = false; args.as_logical = false;
            if (use_vectors)
                args.num_vec_from = &xx_out;
            else
                args.num_pointer_from = xx_out_.get();
        }
        xx_out__ = Rcpp::unwindProtect(SafeRcppVector, (void*)&args);
    }


    return Rcpp::List::create(
        Rcpp::_["ii"] = ii_out__,
        Rcpp::_["jj"] = jj_out__,
        Rcpp::_["xx"] = xx_out__
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
