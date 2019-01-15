# Normalization utility functions
#

#' Utility function to filter out low-expressed genes
#'
#' @params counts
#'
expression_filtering = function(counts, min_counts=20, min_samples=3){
    rows_to_keep = apply(counts, 1, function(x){
            sum(x > min_counts) > min_samples
    })
    counts = counts[rows_to_keep, ]
    return(counts)
}

