library("limma")
library("edgeR")


fit_limma_voom = function(counts, meta, contrasts, use_voom_weights=TRUE){
    y = DGEList(counts=counts)
    y = calcNormFactors(y, method="upperquartile")

    design = model.matrix(~Group + 0, data=meta)

    v = voom(y, design, plot=FALSE)
    v = lmFit(v)

    cleaned_colnames = gsub("Group", "", colnames(design))
    colnames(design) = cleaned_colnames

    allcontrasts = makeContrasts(
        contrasts=contrasts,
        levels=design)

    fit = contrasts.fit(v, allcontrasts)
    fit = topTable(eBayes(fit), number=1000000)

    fit$adj.p.value = p.adjust(fit$P.Value, method="BH")
    return(fit)
}
