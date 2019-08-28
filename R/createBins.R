#' @importFrom GenomicRanges resize width ranges
.tile_variable_regions <- function(gr, bin.width,
                              region.col = 'region_id',
                              region.prefix = 'region_') {
    ## Create variable sized regions, centered on region
    ## that are *each* of width easily divisible by binwidth
    region.widths <- width(ranges(gr))    
    new.widths <- .calc_new_region_widths(region.widths, bin.width)
    gr <- resize(gr, width = new.widths, fix = 'center')

    ## Tile wrt bin.width, annotate, combine into a long GRanges obj
    .tile_annotate_combine(gr, bin.width, region.col, region.prefix)
}

#' @importFrom GenomicRanges resize
.tile_constant_regions <- function(gr, bin.width, region.width,
                                   region.col = 'region_id',
                                   region.prefix = 'region_') {
    ## Center and resize each peak for constant length regions
    new.width <- .calc_new_region_widths(region.width, bin.width)

    ## Produce a message if resizing of region.width occurs
    if (new.width != region.width) {
        message(paste('Automatically resized `region.width` to',
                      new.width,
                      'to be divisible by `bin.width`'))
    }

    ## Center to constant width window around center of peak
    gr <- resize(gr, width = new.width, fix = 'center')

    ## Tile wrt bin.width, annotate, combine into a long GRanges obj
    .tile_annotate_combine(gr, bin.width, region.col, region.prefix)
}






#############################################################
# Internal functions.
#############################################################

#' @importFrom GenomicRanges "mcols<-"
.append_region_id_to_granges <- function(x, row, region.col, region.prefix) {
    ## add a region id column to a GRanges object
    GenomicRanges::mcols(x)[[region.col]] <- paste0(region.prefix, row)
    return(x)
}

#' @importFrom GenomicRanges tile
.tile_annotate_combine <- function(gr, bin.width, region.col, region.prefix) {
    ## Tile each region (peak) by bin width
    tl <- tile(gr, width = bin.width)

    ## Annotate tiles with region id
    tla <- mapply(.append_region_id_to_granges, tl, seq_along(gr),
                  MoreArgs = list(region.prefix = region.prefix,
                                  region.col = region.col))

    ## Combine regions into a single annotated granges and return
    do.call("c", tla)
}

.calc_new_region_widths <- function(region.widths, bin.width) {
    ## make each region width neatly divisible by bin width
    div <- region.widths %% bin.width
    div[div > 0] <- bin.width - div[div > 0]
    new.widths <- region.widths + div
    return(new.widths)
}


#############################################################
# S4 method definitions.
#############################################################

#' @export
setGeneric("createBins", function(x, bin.width = NULL, region.width = NULL, ...) standardGeneric("createBins"))

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
setMethod("createBins", "data.frame", function(x, bin.width = NULL, region.width = NULL, ...) {
    ## Check that prereq columns are present
    req_cols <- c('seqnames', 'start', 'end')
    if (sum(req_cols %in% colnames(x)) != 3) {
        stop("Input must have columns `seqnames`, `start`, and `end`")
    }

    ## Convert to GRanges and run the "default" method
    gr <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)

    createBins(gr, bin.width, region.width, ...)
})


#' @export
setMethod("createBins", "GRanges", function(x, 
                                            bin.width = NULL,
                                            region.width = NULL,
                                            ...) {
    ## Check bin width is provided
    if (is.null(bin.width)) {
        stop('`bin.width` required to be specified')
    }
    
    ## Run the correct method based on the args specified
    if (!is.null(bin.width) & !is.null(region.width)) {
        gr_binned <- .tile_constant_regions(x, bin.width, region.width, ...)
    }
    if (is.null(region.width)) {
        gr_binned <- .tile_variable_regions(x, bin.width, ...)
    }
    
    return(gr_binned)
})


## Testing
## regions_df <- data.frame(seqnames=c('chr1', 'chr2', 'chr1'),
##                   start = c(1, 100, 25), end = c(50, 200, 75))
## regions_gr <- GenomicRanges::makeGRangesFromDataFrame(regions_df)
## bins <- createBins(regions_df, 10)
## createBins(regions_gr, 10)
