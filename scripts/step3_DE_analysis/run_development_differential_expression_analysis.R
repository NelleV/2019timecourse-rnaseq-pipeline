library("limma")
library("edge")
source("utils_edge_contrasts.R")

# FIXME XXX it works, but there's a lot of hacks going on…


###############################################################################
# Options
#

data_dir = "results/varoquaux2019/leaf/"
var_filtering = TRUE
var_cutoff = 0.5
take_log = TRUE

contrasts = c("Control.RT430", "Control.BT642")

###############################################################################
# Load the data
counts = read.delim(file.path(data_dir, "counts_normed.tsv"),
    		    sep="\t",
		    stringsAsFactors=FALSE,
		    check.names=FALSE)
meta = read.table(file.path(data_dir, "meta.tsv"),
		  check.names=FALSE)

if(take_log){
    counts = log2(counts + 1)
}

###############################################################################
# Var filtering
#
# FIXME


###############################################################################
# DE analysis

meta$Time = as.numeric(as.character(meta$Week))
ng_labels = as.factor(
  make.names(meta$Condition:meta$Genotype))
meta$Group = ng_labels


de_analysis = data.frame(row.names=row.names(counts))
for(contrast_formula in contrasts){
    contrast = limma::makeContrasts(
	contrast_formula,
	levels=levels(meta$Group))
    mask = with(meta, Group == contrast_formula)
    model = edgeWithContrasts(as.matrix(counts), meta, contrast, developmental=TRUE,
			      mask=mask)
    contrast_name = gsub(" ", "", contrast_formula, fixed=TRUE)
    contrast_name_pval = paste(contrast_name, "-pval", sep="")
    contrast_name_qval = paste(contrast_name, "-qval", sep="")
    de_analysis[contrast_name_pval] = model$pval
    de_analysis[contrast_name_qval] = model$qval
}

outname = file.path(data_dir, "developmental_de_analysis.tsv")
write.table(
    de_analysis,
    file=outname,
    sep="\t", row.names=TRUE, col.names=TRUE)
