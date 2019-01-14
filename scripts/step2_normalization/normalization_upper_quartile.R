library("config")
library("genefilter")
library("EDASeq")
source("utils_normalization.R")

###############################################################################
# Options
#
# Will filter all_samples that don't have at least `min_counts` counts in
# `min_samples` samples
args = commandArgs(trailingOnly=TRUE)
# For debugging purpose
# data_dir = "results/varoquaux2019/leaf
data_dir = args[1]

config_file = file.path(data_dir, "config.yml")
config = config::get(file=config_file)

filter_expression = unlist(config["filter_expression"] )
min_counts = unlist(config["min_counts"])
min_samples = unlist(config["min_samples"])

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
if(filter_expression){
    counts = expression_filtering(
	counts, meta, min_counts=min_counts, min_samples=min_samples)
}

# Normalize
data_set = EDASeq::newSeqExpressionSet(as.matrix(counts), phenoData=meta)
data_normalized = EDASeq::betweenLaneNormalization(
    data_set, which="upper", offset=TRUE, round=FALSE)
data_normalized = EDASeq::normCounts(data_normalized)


###############################################################################
# Saving the results
write.table(
    data_normalized,
    file=file.path(data_dir, "counts_normed.tsv"),
    sep="\t", row.names=TRUE, col.names=TRUE)
