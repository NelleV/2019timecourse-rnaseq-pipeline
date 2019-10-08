library(moanin)

n_clusters = 2:20
all_labels = list()
for(i in n_clusters){
    filename = paste0("results/stability_", i, ".tsv")
    stability = read.table(filename, sep="\t")
    all_labels[[paste0("B", i)]] = stability
}

png("images/clustering_CDF_consensus.png")
moanin::plot_cdf_consensus(all_labels)
dev.off()

png("images/clustering_AUC_consensus.png")
auc_scores = moanin:::get_auc_similarity_scores(all_labels)
plot(n_clusters, auc_scores, col="black", type="b", pch=16)
dev.off()

png("images/clustering_delta_AUC_consensus.png")
delta_auc = diff(auc_scores)/auc_scores[1:(length(auc_scores)-1)]
plot(3:20, delta_auc, col="black", type="b")

png("images/clustering_AUC_nmi.png")
auc_scores = moanin:::get_auc_similarity_scores(all_labels, method="nmi")
plot(n_clusters, auc_scores, col="black", type="b", pch=16)
dev.off()


png("images/clustering_delta_AUC_nmi.png")
delta_auc = diff(auc_scores)/auc_scores[1:(length(auc_scores)-1)]
plot(3:20, delta_auc, col="black", type="b", pch=16)
dev.off()
