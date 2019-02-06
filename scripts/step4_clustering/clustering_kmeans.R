library(splines)
library(stats)
source("utils_splines.R")
library(devtools)
load_all("../../bin/moanin")

###############################################################################
# Options

n_clusters = 20
df = 6
data_dir = "results/varoquaux2019/leaf/"
filter_genes = TRUE
n_genes_to_keep = 5000
percentage_genes_to_label = 0.5

###############################################################################
# Load the data
counts = read.delim(file.path(data_dir, "counts_normed.tsv"),
    		    sep="\t",
		    stringsAsFactors=FALSE,
		    check.names=FALSE)
meta = read.table(file.path(data_dir, "meta.tsv"),
		  check.names=FALSE)
pvalues = read.table(file.path(data_dir, "timecourse_de_analysis.tsv"),
		     sep="\t",
		     check.names=FALSE)

###############################################################################
# Select "best" genes

if(filter_genes){
    # Filter out q-values for the pvalues table
    col_to_keep = colnames(pvalues)[grepl("-pval", colnames(pvalues))]
    pvalues = pvalues[, col_to_keep]
    fishers_pval = fisher_method(pvalues)
    genes_to_keep = names(sort(fishers_pval)[1:n_genes_to_keep])
    counts = counts[genes_to_keep, ]
}

meta$Time = as.numeric(as.character(meta$Week))
ng_labels = as.factor(
  make.names(meta$Condition:meta$Genotype))
meta$Group = ng_labels

###############################################################################
# Fit the splines

full_model = ~Group:ns(Time, df=df) + Group + 0
X = model.matrix(full_model, data=meta)
y = as.matrix(log2(counts + 1))

fitted_counts = fit_predict_splines(y, X)
rescaled_fitted_counts = rescale_values(fitted_counts, meta)
kmeans_clusters = kmeans(rescaled_fitted_counts, n_clusters, nstart=10)
all_scores = data.frame(row.names=row.names(counts))
for(k in 1:n_clusters){
    scores = score_genes_centroid(y, kmeans_clusters$centers[k,])
    all_scores[paste0("cluster_", k)] = scores
}

# Only assign labels to X% of the genes
max_score = quantile(row_min(all_scores), c(percentage_genes_to_label))
genes_to_not_consider = row_min(all_scores) >= max_score
labels = row_argmin(all_scores)
labels[genes_to_not_consider] = NA
