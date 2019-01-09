library("edge")
library("limma")
source("utils_edge_contrasts.R")

###############################################################################
# Options
#

data_dir = "results/varoquaux2019/leaf/"
var_filtering = TRUE
var_cutoff = 0.5

contrasts = as.vector(
    unlist(read.table(file.path(data_dir, "contrasts"),
		      sep=",", header=FALSE)))

###############################################################################
# Load the data
counts = read.delim(file.path(data_dir, "counts_normed.tsv"),
    		    sep="\t",
		    stringsAsFactors=FALSE,
		    check.names=FALSE)
meta = read.table(file.path(data_dir, "meta.tsv"),
		  check.names=FALSE)

meta$Time = as.numeric(as.character(meta$Week))
ng_labels = as.factor(
  make.names(meta$Condition:meta$Genotype))
meta$Group = ng_labels

de_analysis = data.frame(row.names=row.names(counts))
for(contrast_formula in contrasts){
    contrast = limma::makeContrasts(
	contrast_formula,
	levels=levels(meta$Group))
    model = edgeWithContrasts(as.matrix(counts), meta, contrast)
    contrast_name = gsub(" ", "", contrast_formula, fixed=TRUE)
    contrast_name_pval = paste(contrast_name, "-pval", sep="")
    contrast_name_qval = paste(contrast_name, "-qval", sep="")
    de_analysis[contrast_name_pval] = model$pval
    de_analysis[contrast_name_qval] = model$qval
}

outname = file.path(data_dir, "de_analysis.tsv")
write.table(
    de_analysis,
    file=outname,
    sep="\t", row.names=TRUE, col.names=TRUE)
