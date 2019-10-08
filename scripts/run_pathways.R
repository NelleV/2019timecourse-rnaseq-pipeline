library(devtools)
library(limma)
library(splines)
library(stats)
library(KEGGprofile)
library(kableExtra)
library(viridis)

library(moanin)
library(biomaRt)
df = 6

# Preprocessing & filtering
is_count = FALSE
take_log = FALSE
filter_expression = FALSE

# Clustering
filter_genes = TRUE
n_clusters = 20
percentage_genes_to_label = 0.5
# Set data files and options

data("shoemaker2015")
data = shoemaker2015$data
meta = shoemaker2015$meta

labels = read.table(".results/clustering_labels.txt", sep="\t",
                    check.names=FALSE)

gene_names = row.names(labels)
gene_names = gene_names[!is.na(labels)]
labels = labels[!is.na(labels)]


# FIXME move this to the package
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)

genes = getBM(attributes=c("go_id", "refseq_mrna"),
              values=gene_names,
              filters="refseq_mrna",
              mart=ensembl)

# Create gene to GO id mapping
gene_id_go_mapping = create_go_term_mapping(genes)

clusters = sort(as.vector(unlist(unique(labels))))
all_pathways = data.frame()
for(cluster in clusters){
    genes = gene_names[labels == cluster]

    # convert gene names
    genes = getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
		  filters="refseq_mrna", values=genes,
		  mart=ensembl)["entrezgene_id"]
    genes = as.vector(unlist(genes))
    pathways = KEGGprofile::find_enriched_pathway(
	genes, species="mmu",
	download_latest=TRUE)

    if((dim(all_pathways)[1] == 0) & (dim(pathways$stastic)[1] != 0)){
	print("Initialize")
	pathways$stastic["Cluster"] = cluster
	all_pathways = pathways$stastic
    }else if(dim(pathways$stastic)[1] != 0) {
	print("Merge")
	pathways$stastic["Cluster"] = cluster
	all_pathways = merge(all_pathways, pathways$stastic, all=TRUE)
    }
}


outname = "results/all_pathways.tsv"
write.table(all_pathways, outname, sep="\t")
