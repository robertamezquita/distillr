#' @importFrom csaw regionCounts readParam
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom SummarizedExperiment "rowData<-"
.countBins <- function(x, bins_gr, ..., RCPARAM = readParam(), BPPARAM = SerialParam()) {
    rc <- csaw::regionCounts(x,
                             regions = bins_gr,
                             ...,
                             param = RCPARAM)
    rowData(rc) <- bins_gr
    return(rc)
}    


#############################################################
# Internal functions.
#############################################################


#############################################################
# S4 method definitions.
#############################################################

#' @export
setGeneric("countBins", function(x, bins_gr, ...) standardGeneric("countBins"))

#' @export
setMethod("countBins", "ANY", .countBins)

#' @export
setMethod("countBins", "character", function(x, bins_gr, ...) {
    .countBins(x, bins_gr = bins_gr, ...)
})
