library("config")
library("edge")
library("limma")
library("devtools")
load_all("../../bin/moanin/")

###############################################################################
# Load Options
#

args = commandArgs(trailingOnly=TRUE)
data_dir = args[1]
# For debugging purpose
# data_dir = "results/varoquaux2019/leaf

config_file = file.path(data_dir, "config.yml")
config = config::get(file=config_file)

var_filtering = unlist(config["var_filtering"] )
var_cutoff = unlist(config["var_cutoff"])
take_log = unlist(config["take_log"])

contrasts = as.vector(unlist(
    strsplit(unlist(config["timecourse_contrasts"]), ",")))

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

outname = file.path(data_dir, "timecourse_de_analysis.tsv")
write.table(
    de_analysis,
    file=outname,
    sep="\t", row.names=TRUE, col.names=TRUE)
