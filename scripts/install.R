

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

if(!require(ClusteR)){
    install.packages("ClusteR")
}

if(!require(splines)){
    install.packages("splines")
}

if(!require(NMF)){
    install.packages("NMF")
}

source("https://bioconductor.org/biocLite.R")

if(!require(limma)){
    biocLite("limma")
}

if(!require(edgeR)){
    biocLite("edgeR")
}

if(!require(edge)){
    biocLite("edge")
}

if(!require(EDASeq)){
    biocLite("EDASeq")
}

if(!require(topGO)){
    biocLite("topGO")
}

if(!require(Biobase)){
    biocLite("Biobase")
}

if(!require(MASS)){
    # biocLite("MASS")
}

if(!require("biomaRt")){
    # biocLite("biomaRt")
}

if(!require(KEGGprofile)){
    # biocLite("KEGGprofile")
}

if(!require("NMI")){
    # install.packages("NMI")
}

if(!require(zoo)){
    # install.packages("zoo")
}

if(!require(BiocWorkflowTools)){
    # biocLite("BiocWorkflowTools")
}

