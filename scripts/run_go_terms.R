library(devtools)
library(limma)
library(splines)
library(stats)
library(KEGGprofile)
library(kableExtra)
library(viridis)

library(moanin)
library(biomaRt)
library(topGO)
library(timecoursedata)

df = 6

# Preprocessing & filtering
is_count = FALSE
take_log = FALSE
filter_expression = FALSE

# Differential expression analysis
timecourse_contrasts = c("C-K", "C-M")

# Clustering
filter_genes = TRUE
n_genes_to_keep = 5000
n_clusters = 20
percentage_genes_to_label = 0.5
# Set data files and options
data("shoemaker2015")
data = shoemaker2015$data
meta = shoemaker2015$meta

labels = read.table("results/clustering_labels.txt", sep="\t",
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
ontologies = c("BP")
all_go_terms = data.frame()
for(cluster in clusters){
    for(ontology in ontologies){
        scores = as.numeric(labels == cluster)
        names(scores) = gene_names
        go_terms_enriched = find_enriched_go_terms(
            scores,
            gene_id_go_mapping, ontology=ontology)

        if((dim(all_go_terms)[1] == 0) & (dim(go_terms_enriched)[1] != 0)){
            go_terms_enriched["Cluster"] = cluster
            go_terms_enriched["Ontology"] = ontology
            all_go_terms = go_terms_enriched
        }else if(dim(go_terms_enriched)[1] != 0){
            go_terms_enriched["Cluster"] = cluster
            go_terms_enriched["Ontology"] = ontology
            all_go_terms = merge(all_go_terms, go_terms_enriched, all=TRUE)
        }
    }
}

outname = "results/all_go_terms.tsv"
write.table(all_go_terms, outname, sep="\t")
