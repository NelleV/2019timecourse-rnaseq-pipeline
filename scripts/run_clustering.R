library(moanin)

###############################################################################
# Options

args = commandArgs(trailingOnly=TRUE)
n_clusters = strtoi(args[1])
random_seed = strtoi(args[2])

# Other options
n_genes_to_keep = 5000
filter_genes = TRUE
sample_proportion = 0.8
keep_all_significant = FALSE

outname = paste0(
    "results/stability_", n_clusters,
    "/clustering_", random_seed ,"_labels.tsv")

###############################################################################
# Clustering and everything 
data("shoemaker2015")

data = shoemaker2015$data
meta = shoemaker2015$meta

data = read.table(".results/quantile_normalized.txt") 

# Right now, the package doesn't contain the normalized data
de_analysis = read.table(".results/pvalues.txt", sep="\t", check.names=FALSE)

# Extract pvalues and log_fold change
lfc_col_to_keep = colnames(de_analysis)[
    grepl("-lfc_max", colnames(de_analysis))]
lfc_max = de_analysis[, lfc_col_to_keep]

pval_col_to_keep = colnames(de_analysis)[
    grepl("-pval", colnames(de_analysis))]
pvalues = de_analysis[, pval_col_to_keep]

splines_model = moanin::create_splines_model(meta, degrees_of_freedom=6)

# FIXME refactor this into a function
if(filter_genes){
    # First, filter out any genes that doesn't have a log fold change of at
    # least 2 at at least one of the time points.
    
    genes_to_keep = row.names(lfc_max[rowSums(lfc_max > 2) > 0, ])
    pvalues = pvalues[genes_to_keep, ]

    if(length(genes_to_keep) > n_genes_to_keep){

	# Then rank by fisher's p-value and take max the number of genes of
	# interest
	# Filter out q-values for the pvalues table
	fishers_pval = moanin:::fisher_method(pvalues)
	fishers_qval = stats::p.adjust(fishers_pval)
        if(keep_all_significant){
	   genes_to_keep = names(fishers_qval[fishers_qval < 0.05])
	}else{
	    genes_to_keep = names(sort(fishers_pval)[1:n_genes_to_keep])
	}
    }

    y = as.matrix(data[genes_to_keep, ])
}else{
    y = as.matrix(data)
}

###############################################################################
# Now, randomly subsample.
set.seed(random_seed)
n_genes = dim(y)[1] * sample_proportion
indices = sample(1:dim(y)[1], n_genes, replace=TRUE)

random_subsample_y = y[indices, ]
kmeans_clusters = moanin:::splines_kmeans(
    random_subsample_y, splines_model, n_clusters=n_clusters,
    random_seed=random_seed,
    n_init=20)

# Perform prediction on the whole set of data.
kmeans_clusters = moanin:::splines_kmeans_prediction(
    data, kmeans_clusters)

# Let's save only the labels here.
labels = kmeans_clusters$clusters
if(!dir.exists(dirname(outname))){
    dir.create(dirname(outname))
}
write.table(as.matrix(labels), 
	    outname,
	    sep="\t")
