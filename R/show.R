#' @title Quick Glance at Sparse Objects
#' @description Shows some basic information about a sparse matrix or sparse vector object,
#' without printing a subset of its entries as `Matrix` would do.
#' 
#' Note that this package will by default override the `show` methods of sparse objects, but not
#' the `print` methods - for example, if one defines a variable 'X' containing a sparse matrix, and
#' then types 'X' in the console, that calls the `show` method, but one can still print it by calling
#' 'print(X)'.
#' 
#' In order to restore the `show` method provided by `Matrix`, call `options("MatrixExtra.quick_show" = FALSE)`.
#' @param object A sparse matrix or sparse vector.
#' @param x A sparse vector (same method as in matrix, readded here to avoid naming conflicts).
#' @return The same object that was passed as input, as invisible.
#' @name show
#' @rdname show
#' @examples
#' library(Matrix)
#' library(MatrixExtra)
#' 
#' set.seed(1)
#' X <- Matrix::rsparsematrix(5, 5, .2)
#' set_new_matrix_behavior()
#' show(X)
#' print(X)
#' X
#' 
#' restore_old_matrix_behavior()
#' show(X)
#' print(X)
#' X
NULL

get_sp_type <- function(x) {
    if (inherits(x, "TsparseMatrix")) {
        return("COO")
    } else if (inherits(x, "RsparseMatrix")) {
        return("CSR")
    } else if (inherits(x, "CsparseMatrix")) {
        return("CSC")
    } else {
        throw_internal_error()
    }
}

get_sp_n_entries <- function(x) {
    if (.hasSlot(x, "x")) {
        return(length(x@x))
    } else if (.hasSlot(x, "i")) {
        return(length(x@i))
    } else if (.hasSlot(x, "j")) {
        return(length(x@j))
    } else {
        throw_internal_error()
    }
}

## avoid integer overflow when calculating matrix dimensions
get_pct_full <- function(x, nnz) {
    temp <- as.numeric(nnz) / as.numeric(nrow(x))
    temp <- temp / ncol(x)
    return(temp * 100)
}

autoadj_pct <- function(x) {
    if (is.na(x) || !is.finite(x)) return(as.character(x))
    if (x >= 1) return(sprintf("%.2f%%", x))

    places <- -floor(log10(x))
    if (is.na(places) || !is.finite(places)) {
        places <- 2
    }
    places <- max(places, 2)
    places <- min(places, 7)
    fmt <- sprintf("%%.%df%%%%", places)
    return(sprintf(fmt, x))
}

prettyNum_custom <- function(x) {
    if (inherits(x, "integer")) {
        x <- sprintf("%d", x)
    } else {
        x <- sprintf("%.0f", x)
    }
    x <- prettyNum(x, big.mark=",")
    return(x)
}

glance_sp_matrix <- function(object) {
    if (!getOption("MatrixExtra.quick_show", default=FALSE)) {
        Matrix::printSpMatrix2(object)
        return(invisible(object))
    }

    if (length(class(object)) == 1L) {
        cat(sprintf("Sparse %s matrix (class '%s')\n", get_sp_type(object), class(object)))
    } else {
        classes <- paste(class(object), sep=", ")
        cat(sprintf("Sparse %s matrix (classes %s)\n", get_sp_type(object), class(object)))
    }
    cat(sprintf("Dimensions: %s x %s\n",
        prettyNum_custom(nrow(object)),
        prettyNum_custom(ncol(object))))
    if (nrow(object) && ncol(object)) {
        nnz <- get_sp_n_entries(object)
        cat(sprintf("(%s entries, %s full)\n", prettyNum_custom(nnz),
            autoadj_pct(get_pct_full(object, nnz))))
    }
    return(invisible(object))
}

#' @rdname show
#' @export
setMethod("show", signature(object="sparseMatrix"), glance_sp_matrix)


### Taken from Matrix:
### https://github.com/cran/Matrix/blob/2e9cdfc101b37e0dfe0ca046afdacb1d40b720c4/R/sparseVector.R
###
### These methods are not exported so they were copy-pasted here in order to be able to
### switch between overriding them and not overriding them
prSpVector <- function(x, digits = getOption("digits"),
            maxp = getOption("max.print"), zero.print = ".")
{
    cld <- getClassDef(class(x))
    stopifnot(extends(cld, "sparseVector"), maxp >= 1)
    if(is.logical(zero.print))
    zero.print <- if(zero.print) "0" else " "
##     kind <- .M.kindC(cld)
##     has.x <- kind != "n"
    n <- x@length
    if(n > 0) {
        if(n > maxp) { # n > maxp =: nn : will cut length of what we'll display :
            x <- head(x, maxp)
            n <- maxp
        }
        xi <- x@i
        is.n <- extends(cld, "nsparseVector")
        logi <- is.n || extends(cld, "lsparseVector")
        cx <- if(logi) rep.int("N", n) else character(n)
        cx[if(length(xi)) -xi else TRUE] <- zero.print
        cx[ xi] <- {
        if(is.n) "|" else if(logi) c(":","|")[x@x + 1L] else
        ## numeric (or --not yet-- complex): 'has.x' in any cases
        format(x@x, digits = digits)
        }
        ## right = TRUE : cheap attempt to get better "." alignment
        print(cx, quote = FALSE, right = TRUE, max = maxp)
    }
    invisible(x) # TODO? in case of n > maxp, "should" return original x
}

old_show_spvec <- function(object) {
    n <- object@length
    cl <- class(object)
    cat(sprintf('sparse vector (nnz/length = %d/%.0f) of class "%s"\n',
        length(object@i), as.double(n), cl))
    maxp <- max(1, getOption("max.print"))
    if(n <= maxp) {
        prSpVector(object, maxp = maxp)
    } else { # n > maxp : will cut length of what we'll display :
        ## cannot easily show head(.) & tail(.) because of "[1] .." printing of tail
        prSpVector(head(object, maxp), maxp = maxp)
        cat(" ............................",
            "\n ........suppressing ", n - maxp,
            " entries in show(); maybe adjust 'options(max.print= *)'",
            "\n ............................\n\n", sep='')
   }
   invisible(object)
}


glance_sp_vector <- function(object) {
    if (!getOption("MatrixExtra.quick_show", default=FALSE)) {
        old_show_spvec(object)
        return(invisible(object))
    }

    if (length(class(object)) == 1L) {
        cat(sprintf("Sparse vector (class '%s')\n", class(object)))
    } else {
        classes <- paste(class(object), sep=", ")
        cat(sprintf("Sparse vector (classes: %s)\n", classes))
    }
    cat(sprintf("Length: %s\n", prettyNum_custom(object@length)))
    if (length(object)) {
        nnz <- get_sp_n_entries(object)
        cat(sprintf("(%s entries, %s full)\n",
            prettyNum_custom(nnz),
            autoadj_pct(100.0 * (nnz / length(object)))))
    }
    return(invisible(object))
}

#' @rdname show
#' @export
setMethod("show", signature(object="sparseVector"), glance_sp_vector)

#' @rdname show
#' @export
setMethod("print", signature(x="sparseVector"), function(x) old_show_spvec(x))
