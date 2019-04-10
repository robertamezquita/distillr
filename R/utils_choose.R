## lifted from scater: R/plot_utils.R

#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom SingleCellExperiment int_elementMetadata int_colData
#' @importFrom BiocGenerics rownames colnames
.choose_values <- function(x, by, mode=c("column", "row"), search=c("any", "metadata", "exprs"),
    exprs_values = "logcounts", coerce_factor = FALSE, level_limit = Inf, discard_solo = FALSE) 
# This function looks through the SCE data and returns the
# values to be utilized. Either 'by' itself, or a column of colData,
# or a column of rowData (or the expression values of a feature)
{
    vals <- NULL
    if (is.character(by) && length(by) > 0) { 
        mode <- match.arg(mode)
        search <- match.arg(search)
        internal_only <- FALSE

        # Determining what to check, based on input 'by'.
        if (length(by)>1) {
            if (search=="exprs") {
                stop("character vector of length > 1 not allowed for search='exprs'")
            }
            search <- "metadata"

            if (is.na(by[1])) { 
                internal_only <- TRUE
                by <- by[-1]
            }
        } else if (search=="any") { 
            cur_name <- names(by)
            if (!is.null(cur_name) && !is.na(cur_name)) { 
                if (cur_name=="metadata" || cur_name =="exprs") {
                    search <- cur_name
                } 
            }
        }
        names(by) <- NULL
           
        # Checking the metadata; note the loop to account for nesting.
        if (search=="any" || search=="metadata") {
            if (!internal_only) { 
                meta_data <- if (mode=="column") colData(x) else rowData(x)
            } else {
                meta_data <- if (mode=="column") int_colData(x) else int_elementMetadata(x)
            }
            for (field in by) {
                if (!field %in% colnames(meta_data)) {
                    break
                }
                vals <- meta_data[[field]]
                meta_data <- vals
            }
            by <- paste(by, collapse=":") # collapsing to a single string for output.
        }

        # Metadata takes priority, so we don't bother searching if 'vals' is non-NULL.
        if ((search=="any" || search=="exprs") && is.null(vals)) { 
            exprs <- assay(x, i = exprs_values, withDimNames=FALSE)
            if (mode=="column") {
                m <- match(by, rownames(x)) # coloring columns, so we take the row values.
                if (!is.na(m)) {
                    vals <- exprs[m,] 
                }
            } else if (mode=="row") {
                m <- match(by, colnames(x)) 
                if (!is.na(m)) {
                    vals <- exprs[,m] 
                }
            }
        }

        if (is.null(vals)) {
            stop(sprintf("cannot find '%s' in %s fields", by, search))
        }

    } else if (is.data.frame(by)) {
        if (ncol(by) != 1L) {
            stop("input data frame should only have one column")
        } else {
            if (mode=="column" && nrow(by) != ncol(x)) {
                stop("number of rows of input data frame should be equal to 'ncol(object)'")
            }
            if (mode=="row" && nrow(by) != nrow(x)) {
                stop("number of rows of input data frame should be equal to 'nrow(object)'")
            }
        }

        ## Allow arbitrary values to be specified.
        vals <- by[,1]
        by <- colnames(by)

    } else if (!is.null(by)) {
        # We have to allow by=NULL to pass through smoothly.
        stop("invalid value of 'by' supplied")
    }

    # Checking the level limit.
    if (coerce_factor && !is.null(vals)) {
        vals <- factor(vals)
        if (level_limit < nlevels(vals)) {
            stop(sprintf("number of unique levels exceeds %i", level_limit))
        }
    }

    # If only one level for the variable, set to NULL.
    if (length(unique(vals))<=1L && discard_solo) { 
        by <- NULL
        vals <- NULL
    }
    return(list(name = by, val = vals))
}
