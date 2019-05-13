


folders = list.files("results/")
folders = folders[grepl("stability", folders) & !grepl("tsv", folders)]

for(folder in folders){
    print(folder)
    filenames = list.files(paste0("results/", folder))
    all_labels = NULL
    for(filename in filenames){
	labels = unique(
	    read.table(paste0("results/", folder, "/", filename),
		       row.names=NULL))
	colnames(labels) = c("Gene", "Label")
	row.names(labels) = labels$Gene
	labels = subset(labels, select=c("Label"))
	all_labels = merge(all_labels, labels, by="row.names", all=TRUE)
	row.names(all_labels) = all_labels$Row.names
	all_labels = subset(all_labels, select=-c(Row.names))
    }
    colnames(all_labels) = lapply(
	1:dim(all_labels)[2],
	function(x){paste0("C", x)})
    outname = paste0("results/", folder, ".tsv")
    write.table(all_labels, outname, sep="\t")
}
