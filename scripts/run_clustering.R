library(moanin)

###############################################################################
# Options

args = commandArgs(trailingOnly=TRUE)
n_clusters = strtoi(args[1])
random_seed = strtoi(args[2])

# Other options
filter_genes = TRUE
sample_proportion = 1
keep_all_significant = FALSE

outname = paste0(
    "results/stability_", n_clusters,
    "/clustering_", random_seed ,"_labels.tsv")

###############################################################################
# Clustering and everything 
data("shoemaker2015")

data = shoemaker2015$meta
meta = shoemaker2015$meta

# Right now, the package doesn't contain the normalized data
de_analysis = read.table(".results/pvalues.txt", sep="\t", check.names=FALSE)

# Extract pvalues and log_fold change
lfc_col_to_keep = colnames(de_analysis)[
    grepl("-lfc_max", colnames(de_analysis))]
lfc_max = de_analysis[, lfc_col_to_keep]

pval_col_to_keep = colnames(de_analysis)[
    grepl("-pval", colnames(de_analysis))]
pvalues = de_analysis[, pval_col_to_keep]

splines_model = moanin::create_moanin_model(meta, degrees_of_freedom=6)
fishers_pval = moanin:::pvalues_fisher_method(pvalues)
fishers_qval = stats::p.adjust(fishers_pval)

genes_to_keep = row.names(
    lfc_max[(rowSums(lfc_max > 2) > 0) & (fishers_qval < 0.05), ])
y = as.matrix(data[genes_to_keep, ])


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
kmeans_clusters = moanin::: splines_kmeans_predict(
    y, kmeans_clusters)

# Let's save only the labels here.
labels = kmeans_clusters$clusters
if(!dir.exists(dirname(outname))){
    dir.create(dirname(outname))
}
write.table(as.matrix(labels), 
	    outname,
	    sep="\t")
