library(reshape2)


#' Compute consensus matrix from labels
#'
#' @param labels matrix (n_genes, n_clusters)
#'	Matrix containing the set of labels
#' @param scale	boolean, optional, default: TRUE
#'	Whether to rescale the matrix
#' @export
consensus_matrix = function(labels, scale=TRUE){
    melted_labels = reshape2::melt(as.matrix(labels))
    colnames(melted_labels) = c("Gene", "Clustering", "Label")
    melted_labels = subset(melted_labels, select=c("Gene", "Label"))
    melted_labels = melted_labels[!is.na(melted_labels["Label"]),]
    consensus = tcrossprod(table(melted_labels)) 
    if(scale){
        diag_elements = diag(consensus)
	consensus = t(t(consensus / diag_elements) / diag_elements)
	diag(consensus) = 0
    }
    return(consensus)
}
