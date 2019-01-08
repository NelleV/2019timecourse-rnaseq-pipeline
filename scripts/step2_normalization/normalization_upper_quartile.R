library("genefilter")
library("EDASeq")

###############################################################################
# Options
#
# Will filter all_samples that don't have at least `min_counts` counts in
# `min_samples` samples
data_dir = "results/varoquaux2019/leaf/"
expression_filtering = TRUE
min_counts = 20
min_samples = 3

###############################################################################
# Load the data
counts = read.delim(file.path(data_dir, "data.tsv"),
    		    sep="\t",
		    stringsAsFactors=FALSE,
		    check.names=FALSE)
meta = read.table(file.path(data_dir, "meta.tsv"),
		  stringsAsFactors=FALSE,
		  check.names=FALSE)

###############################################################################
# Filter very low expressed genes before normalization
expression_filtering = function(counts, meta, min_counts=20, min_samples=3){
    rows_to_keep = apply(counts, 1, function(x){
            sum(x > min_counts) > min_samples
    })
    counts = counts[rows_to_keep, ]
    return(counts)
}
counts = expression_filtering(counts, meta, min_counts=min_counts, min_samples=min_samples)

# Normalize
data_set = EDASeq::newSeqExpressionSet(as.matrix(counts), phenoData=meta)
data_normalized = EDASeq::betweenLaneNormalization(data_set, which="upper", offset=TRUE, round=FALSE)
data_normalized = EDASeq::normCounts(data_normalized)


###############################################################################
# Saving the results
write.table(
    data.frame("gene_id"=row.names(data_normalized),
	       data_normalized, check.names=FALSE),
    file=file.path(data_dir, "counts_normed.txt"),
    sep="\t", row.names=FALSE, col.names=TRUE)
