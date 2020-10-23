

# This is a bit annoying, but I'm hoping to save time on travis-ci
if(!require("devtools")){
    install.packages("devtools")
}

if(!require(ggfortify)){
    install.packages("ggfortify")
}

if(!require(kableExtra)){
    install.packages("kableExtra")
}

if(!require(testthat)){
    install.packages("testthat")
}

if(!require(viridis)){
    install.packages("viridis")
}

if(!require(knitr)){
    install.packages("knitr")
}

if(!require(rmarkdown)){
    install.packages("rmarkdown")
}

if(!require(roxygen2)){
    install.packages("roxygen2")
}

if(!require(ClusterR)){
    install.packages("ClusterR")
}

if(!require(splines)){
    install.packages("splines")
}

if(!require(NMF)){
    install.packages("NMF")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("topGO")

if(!require(limma)){
    BiocManager::install("limma")
}

if(!require(edgeR)){
    BiocManager::install("edgeR")
}

if(!require(edge)){
    BiocManager::install("edge")
}

if(!require(EDASeq)){
    BiocManager::install("EDASeq")
}

if(!require(topGO)){
    BiocManager::install("topGO")
}

if(!require(Biobase)){
    BiocManager::install("Biobase")
}

if(!require(MASS)){
    BiocManager::install("MASS")
}

if(!require("biomaRt")){
    BiocManager::install("biomaRt")
}

if(!require(KEGGprofile)){
    BiocManager::install("KEGGprofile")
}

if(!require("NMI")){
    install.packages("NMI")
}

if(!require(zoo)){
    install.packages("zoo")
}

if(!require(BiocWorkflowTools)){
    BiocManager::install("BiocWorkflowTools")
}

if(!require("pander")){
    install.packages('pander')
}

if(!require(rticles)){
    install.packages("rticles")
}

if(!require(matrixStats)){
    install.packages("matrixStats")
}
