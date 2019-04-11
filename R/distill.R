#' @importFrom BiocParallel bplapply SerialParam
.distill <- function(x,
                     compare, groups,
                     regions, 
                     band = 5,
                     quantile = 0.9,
                     adjust_variance = TRUE, 
                     subset.row = NULL, subset.col = NULL,
                     BPPARAM = SerialParam())
{
    ## `groups` must correspond to columns and describe their respective groups
    ##    - e.g. a vector of group labels
    ## `compare` gives two groups to specify comparison to be made
    ##    - e.g. a length 2 character/factor
    ## `regions` must correspond to rows and describe their respective grouping
    ##    - e.g. a vector of region labels
    ## `band`/`quantile` are for the statistical method

    ## TODO: use compare arg to order the comparison (should change sign)
    
    ## Subsetting --------------------------------------------------------------
    ## Check and perform subsetting (returns numeric indices)
    subset.row <- .subset_to_index(subset.row, x, byrow = TRUE)
    subset.col <- .subset_to_index(subset.col, x, byrow = FALSE)
    x <- x[subset.row, subset.col]
    
    ## Subset the groups/regions variables to match data in `x` by cols/rows
    regions <- regions[subset.row]    
    groups <- groups[subset.col]

    ## Check input of `compare` arg
    ## Check compare (after group subsetting)
    if (length(compare) != 2) {
        stop("`compare` must be of length 2")
    }
    if (sum(compare %in% groups) != 2) {
        stop("Specified groups in `compare` not found in `groups`")
    }

    ## Subset down to what is being compared
    xc <- x[, groups %in% compare]
    groups_xc <- as.character(groups)[groups %in% compare]


    ## Calc combinations to compare and make result names ----------------------
    cn <- combn(1:length(groups_xc), 2) # cols = uniq comb, row = g1 + g2

    ## Create names for test statistics output for each comparison
    if (!is.null(colnames(xc))) {
        grabbed_colnames <- apply(cn, 2, function(x) { colnames(xc)[x] })
        ts_names <- apply(grabbed_colnames, 2, paste0, collapse = '..')
    } else {
        grabbed_cols <- apply(cn, 2, paste0, collapse = '..')
        gn <- combn(groups_xc, 2)           # comparison made        
        gnp <- apply(gn, 2, paste0, collapse = '..') # comparison made, labelized        
        ts_names <- paste0(grabbed_cols, '_', gnp)
    }

    ## Calculate test statistic per comparison ---------------------------------
    ## Reformat data into list: samples then regions
    xl <- apply(xc, 2, split, regions)

    ## Split by comparison
    ts_df <- list()
    for (i in 1:ncol(cn)) {
        s1 <- xl[[cn[1, i]]]
        s2 <- xl[[cn[2, i]]]

        ## Split by region per comparison
        ## Estimate the variance in data
        sigma <- bpmapply(.core_estimate_variance, s1, s2,
                          MoreArgs = list(band = band),
                          BPPARAM = BPPARAM)

        ## Adjust the variance using an empirical bayesian approach
        ## Any variances that are equal to zero are set to min nonzero value
        sigma[sigma == 0] <- min(sigma[sigma > 0])

        ## Adds to the variances +quantile of calculate variances to "play it safe"
        ## against low variance regions
        if (adjust_variance == TRUE) {
            sigma <- sigma + quantile(sigma, quantile)
        }
        sigma_l <- as.list(sigma) # convert to list for mapply func
        
        ## Test statistics portion
        ## Calculate the squared kernel estimator for T_lambda
        ## we get this under the data not the null hypothesis
        ts_l <- bpmapply(.core_signed_ts, s1, s2, sigma_l,
                         MoreArgs = list(band = band), BPPARAM = BPPARAM)

        ## Convert to data frame ordered by region appearance; assign into list
        tstat <- c(unlist(ts_l[1, ]), use.names = names(ts_l[1, ]))
        tsign <- c(unlist(ts_l[2, ]), use.names = names(ts_l[2, ]))
        tmean <- c(unlist(ts_l[3, ]), use.names = names(ts_l[3, ]))

        ts_df[[i]] <- data.frame(Tstat = tstat, Tsign = tsign, Tmean = tmean)[unique(regions), ]
    }
    names(ts_df) <- ts_names

    ## Return the list of data frames of test statistics per comparison performed
    return(ts_df)
}


#############################################################
# Internal functions.
#############################################################

.create_hat_matrix <- function(n, band) {
    ## Create Nadaraya-Watson Hat matrix
    ## - square, sym matrix = individual kernels corresponding
    ##   to normal dist centered at each bin, where increasing
    ##   i moves along the length of the region
    ind <- 1:n / n        # bin indicator var
    hwidth <- band / n    # corollary for proportion of data observed per bin
    
    Idt <- diag(rep(1, n))
    Snw <- apply(Idt, 2, function(y) {
        ksmooth(ind, y, 
                kernel = 'normal', bandwidth = hwidth,
                x.points = ind)$y
    })
    return(Snw)
}

.core_estimate_variance <- function(r1, r2, band) {
    diff <- r1 - r2       # difference vector
    n <- length(diff)     # number of bins observed
    hwidth <- band / n    # corollary for proportion of data observed per bin
    ind <- 1:n / n        # bin indicator var

    ## Create the Nadaraya-Watson hat matrix
    Snw <- .create_hat_matrix(n, band)
    
    ## Get kernel smoothed version of difference vector
    diff_ks <- ksmooth(ind, diff, kernel = 'normal', bandwidth = hwidth)$y


    ## Calculate degrees of freedom
    ## ad_df <- 2 * sum(diag(Snw)) - sum(diag(Snw %*% t(Snw))) ## closed form
    ad_df <- sum(diag(Snw)) # degrees of freedom

    ## Calculate sigma (RSS) from the difference between smoothed v raw version
    sigma <- sqrt(sum((diff_ks - diff)^2) / (length(diff_ks) - ad_df))

    return(sigma)
}

    
.core_signed_ts <- function(r1, r2, sd, band) {
    diff <- r1 - r2       # difference vector
    n <- length(diff)     # number of bins observed
    hwidth <- band / n    # corollary for proportion of data observed per bin
    ind <- 1:n / n        # bin indicator var

    ## data now adjusted by new sd
    diff_adj <- diff / sd

    ## Get kernel smoothed version of diff vector using adj. sd
    diff_adj_ks <- ksmooth(ind, diff_adj, kernel = 'normal', bandwidth = hwidth)$y

    ## Get mean and sign of overall test statistic
    ts_mean <- mean(diff_adj_ks^2)
    sign <- sign(sum(diff_adj_ks))

    ## Get final test stat under the null with the Wilson-Hilferty transform
    Snw <- .create_hat_matrix(n, band)
    Amax <- Snw %*% Snw
    eigenvalue <- as.numeric(eigen(Amax)$values)
    dof <- sum(eigenvalue)^2 / sum(eigenvalue^2)
    delta <- sum(eigenvalue) / (n * dof)
    
    ts_kn <- ((ts_mean / (delta * dof))^(1/3) - (1 - 2 / (9 * dof))) / (sqrt(2 / (9 * dof)))
    
    ## Combine results and return ----------------------------------------------
    return(list(Tstat = ts_kn, Tsign = sign, Tmean = ts_mean))
}


#############################################################
# S4 method definitions.
#############################################################

#' @export
setGeneric("distill", function(x, compare, groups, regions, band, quantile, ...) standardGeneric("distill"))

#' @export
setMethod("distill", "ANY", .distill)


#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment "rowData<-" rowData
#' @importFrom S4Vectors "metadata<-"
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom S4Vectors cbind
setMethod("distill", "SummarizedExperiment", function(x,
                                                      compare, group_by,
                                                      region_by,
                                                      band, quantile,
                                                      ...,
                                                      assay.type = "counts",
                                                      res.out = FALSE)
{
    ## `group_by` should refer to a column
    ## ... should *not* include:
    ##   - groups (calculated from `group_by` for SCEs)
    ##   - regions (calculated from `region_by` for SCEs from rowData())
    ## ... should *def* include:
    ##   - band/quantile for stat method
    ## ... should *possibly* include:
    ##   - subset.row/subset.col
    ##   - BPPARAM()
    ## `region_by` should refer to a column in metadata that groups regions

    ## Extract grouping variable values
    groups <- .choose_values(x, group_by, mode = 'column', search = 'metadata')$val
    regions <- .choose_values(x, region_by, mode = 'row', search = 'metadata')$val

    ## Coerce grouping vars to character
    groups <- as.character(groups)
    regions <- as.character(regions)
    
    ## Run main method
    ds <- .distill(assay(x, i = assay.type),
                   compare = compare, groups = groups,
                   regions = regions,
                   band, quantile,
                   ...)

    ## Shortcircuit results wrangling - output results directly
    if (res.out == TRUE) {
        return(ds)
    }

    ## TEMP: return with results appended to metadata
    ## TODO: return with results added to rowData as a list (compact = TRUE)
    ##       - ideally would be a nested DataFrame
    ## TODO: (?) engineer a new struct, rowGroupData, to organize rowData <-> groups relations
    out_name <- paste0('distillr_',
                       paste0(compare, collapse = '..'),
                       '_band-', band,
                       '_quantile-', quantile)
    
    metadata(x)[[out_name]] <- ds

    return(x)
})





        ## SCRAP ====================

        
    ##         ## TEMP RETURN DUMMY OBJECT ------------------------------------------------
    ## ## Calc how many region (row) groups there are
    ## row_groups <- length(unique(regions))

    ## ## TEMP: use the BPPARAM argr
    ## tmp <- bplapply(1:10, sum, BPPARAM = BPPARAM)
    
    ## ## Return dummy object
    ## out <- list(
    ##     ## per row stats (for rowData/rowRanges)
    ##     per_row_stats = DataFrame(diff = rep(0, nrow(x)),
    ##                               hold = rep(1, nrow(x))),
    ##     ## row stats - nested/grouped by region (for metadata)
    ##     per_rowgroup_stats = DataFrame(stat = rep(0, row_groups),
    ##                                    pval = rep(0, row_groups)),
    ##     ## parametrization
    ##     method_params = list(compare = compare,
    ##                          band = band,
    ##                          quantile = quantile,
    ##                          subset.row = subset.row,
    ##                          subset.col = subset.col)
    ## )
