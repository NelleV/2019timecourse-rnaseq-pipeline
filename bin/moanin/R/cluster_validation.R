library(reshape2)


consensus_matrix = function(labels){
    melted_labels = reshape2::melt(labels)
    colnames(melted_labels) = c("Gene", "Clustering", "Label")
    melted_labels = subset(melted_labels, select=c("Gene", "Label"))
    return(tcrossprod(table(melted_labels)))
}
