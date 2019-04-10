#' @importFrom BiocGenerics cbind
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
.convert_to_nested_DataFrame <- function(existing, set_list, stat_list, exprs_values = "counts") {
    n_values <- length(stat_list[[1]][[1]]) # There should be at least one statistic.
    output <- .create_outer_DataFrame(set_list, n_values)

    sub_output <- new("DataFrame", nrows=n_values)
    for (x in names(stat_list)) { 
        current <- stat_list[[x]] # need to do it via "[[<-" to store DataFrames as columns.
        current <- .cbind_overwrite_DataFrames(existing[[x]], current)
        sub_output[[x]] <- current
    }
    
    output <- cbind(output, sub_output)
    .cbind_overwrite_DataFrames(existing, output)
}

#' @importFrom BiocGenerics colnames<- cbind
.convert_to_full_DataFrame <- function(existing, set_list, stat_list, trim.fun=identity) {
    n_values <- length(stat_list[[1]][[1]]) # There should be at least one statistic.
    output <- .create_outer_DataFrame(set_list, n_values)

    collected <- stat_list
    for (x in names(stat_list)) { 
        current <- stat_list[[x]]

        # For consistency with old output.
        if (x!="all") { 
            colnames(current) <- sprintf("%s_%s", colnames(current), trim.fun(x))
        }

        collected[[x]] <- current
    }

    combined <- do.call(cbind, c(list(output), unname(collected)))
    .cbind_overwrite_DataFrames(existing, combined)
}

#' @importFrom S4Vectors DataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
.create_outer_DataFrame <- function(set_list, n_values) {
    if (length(set_list)) { 
        output <- do.call(DataFrame, set_list)
        colnames(output) <- sprintf("is_%s", colnames(output)) 
    } else {
        output <- new("DataFrame", nrows=n_values)
    }
    return(output)
}

#' @importFrom BiocGenerics colnames cbind
.cbind_overwrite_DataFrames <- function(existing, updated) {
    if (is.null(existing)) { 
        return(updated)
    }
    existing <- existing[, !(colnames(existing) %in% colnames(updated)), drop = FALSE]
    cbind(existing, updated)   
} 
