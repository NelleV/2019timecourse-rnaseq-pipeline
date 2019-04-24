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
    melted_labels$Clustering_lab = melted_labels$Clustering:as.factor(melted_labels$Label)
    w = reshape2::dcast(melted_labels, Gene~Clustering_lab)
    x = as.matrix(w[,-1])
    x[is.na(x)] = 0
    x = apply(x, 2,  function(x) as.numeric(x > 0))
    consensus = tcrossprod(x) 
    if(scale){
        diag_elements = diag(consensus)
	consensus = t(t(consensus / diag_elements) / diag_elements)
	diag(consensus) = 0
    }
    return(consensus)
}
