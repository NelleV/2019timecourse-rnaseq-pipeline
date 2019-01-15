library("limma")
source("utils_limma_voom.R")


###############################################################################
# Options
#
#args = commandArgs(trailingOnly=TRUE)
#data_dir = args[1]
# For debugging purpose
data_dir = "results/varoquaux2019/leaf"

config_file = file.path(data_dir, "config.yml")
config = config::get(file=config_file)

var_filtering = unlist(config["var_filtering"] )
var_cutoff = unlist(config["var_cutoff"])
take_log = unlist(config["take_log"])
contrasts = as.vector(unlist(strsplit(unlist(config["weekly_contrasts"]), ",")))


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

de_analysis = fit_limma_voom(counts, meta, contrasts)

outname = file.path(data_dir, "weekly_de_analysis.tsv")
write.table(
    de_analysis,
    file=outname,
    sep="\t", row.names=TRUE, col.names=TRUE)
