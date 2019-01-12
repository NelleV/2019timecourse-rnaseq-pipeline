library("limma")
source("utils_limma_voom.R")


###############################################################################
# Options
#

data_dir = "results/varoquaux2019/leaf/"
var_filtering = TRUE
var_cutoff = 0.5
take_log = TRUE

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

for(w in unique(meta$Week)){
    for(contrast_formula in contrasts){
        contrast_name = gsub(" ", "", contrast_formula, fixed=TRUE)
	contrast_name_pval = paste(contrast_name, "-timepoint", w, "-pval", sep="")
	contrast_name_qval = paste(contrast_name, "-timepoint", w, "-qval", sep="")
	contrast_name_lfc = paste(contrast_name, "-timepoint", w, "-lfc", sep="")

        mask = with(meta, Week == w)
	c = counts[, mask]
	m = meta[mask, ]
	tryCatch(
	    {
		model = fit_limma_voom(c, m, contrast_formula)
	    },
	    error=function(err) { model = NULL},
	    finally={
		de_analysis[contrast_name_pval] = model$P.Value
		de_analysis[contrast_name_qval] = model$adj.P.Val
		de_analysis[contrast_name_lfc] = model$logFC
	    }
	)
    }
}

outname = file.path(data_dir, "weekly_de_analysis.tsv")
write.table(
    de_analysis,
    file=outname,
    sep="\t", row.names=TRUE, col.names=TRUE)
