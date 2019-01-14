
###############################################################################
# Options
#

data_dir = "data/varoquaux2019/root/"
results_dir = "results/varoquaux2019/root/"

###############################################################################
# Load the data
counts = read.delim(file.path(data_dir, "data.tsv"),
    		    sep="\t",
		    stringsAsFactors=FALSE,
		    check.names=FALSE)
meta = read.table(file.path(data_dir, "meta.tsv"),
		  check.names=FALSE)

###############################################################################
# The pipeline requires:
#   - a Group column in the metadata file, containing the "group" used to do
#     the time-course differential expression analysis
#   - a Time colmun, containing the time variable (in the case of EPICON, the
#     week
meta$Time = as.numeric(as.character(meta$Week))
ng_labels = as.factor(
  make.names(meta$Condition:meta$Genotype))
meta$Group = ng_labels


###############################################################################
# Now, just save in the folder
outname = file.path(results_dir, "data.tsv")
write.table(counts, file=outname, sep="\t", row.names=TRUE, col.names=TRUE)

outname = file.path(results_dir, "meta.tsv")
write.table(meta, file=outname, sep="\t", row.names=TRUE, col.names=TRUE)

file.copy(file.path(data_dir, "contrasts"),
	  file.path(results_dir, "contrasts"))

file.copy(file.path(data_dir, "config.yml"),
	  file.path(results_dir, "config.yml"))
